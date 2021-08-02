__author__  =  'Jon K Peck'
__version__ =  '1.0'
version = __version__

# history
# 07-31-2021 initial version

# The STATS NTILE ANALYSIS extension command

import spss, spssaux, spssdata
from extension import Template, Syntax, processcmd
import bisect, random

# debugging
        # makes debug apply only to the current thread
try:
    import wingdbstub
    import threading
    wingdbstub.Ensure()
    wingdbstub.debugger.SetDebugThreads({threading.get_ident(): 1})
except:
    pass


# main routine
def dontile(predprob, depvar, ntiles=10, predval=1.,
        gainchart=False, responsechart=False, liftchart=False):
    """See Run for argument definitions"""
    
    if spss.GetSplitVariableNames():
        raise ValueError(_("The STATS NTILE ANALYSIS procedure does not support SPLIT FILES."))
    
    mainds = spss.ActiveDataset()
    anycharts = any([gainchart, responsechart, liftchart])
    if mainds == "*" and anycharts:
        raise ValueError(_("The active dataset must have a name in order to use this procedure if charting"))
    
    # the target value is a string here but needs to be numeric if the depvar is numeric
    vardict = spssaux.VariableDict(depvar)
    isstrvar = vardict[depvar].VariableType > 0
    vallbls = vardict[depvar].ValueLabelsTyped
    if not isstrvar:
        predval = float(predval)
    targetlabel = vallbls.get(predval, None)
    if targetlabel is None:
        targetlabel = str(predval)
    else:
        targetlabel = f"{predval} ({targetlabel})"
    # get the ntiles of the predicted probabilities
    ntilelist, ntileupperbnds, nbnds = dontiles(ntiles, predprob)
    
    cases = {}
    targetcount = {}
    minp = {}
    maxp = {}
    # a table cell might wind up empty so initial all possible cells to 0.
    for i in range(len(ntilelist)):
        cases[i] = 0.
        targetcount[i] = 0.
        minp[i] = None
        maxp[i] = None
        
    casevars = [predprob, depvar]
    wtvar = spss.GetWeightVar()
    if wtvar:
        casevars.append(wtvar)
    wt = 1.
    
    # pass the data and accumulate interval statistics
    curs = spssdata.Spssdata(casevars, names=False, omitmissing=True)
    try:
        for case in curs:
            if wtvar:
                wt = case[2]
            if wt <= 0.:   #zero-weighted and negative-weighted cases are not seen by procedures
                continue
            
            # The xml table output from FREQ only has 14 decimals, but the variable values
            # can have a few more, which could exceed the last actual upperbound
            # so we gently round the predicted probabilities a bit
            tallyval = round(case[0], 13)
            loc = bisect.bisect_left(ntileupperbnds, tallyval)  # find the interval for the case
            if loc >= nbnds:
                raise ValueError(_(f"A probability in the data exceeds 1: {tallyval}"))
            cases[loc] += wt
            if minp[loc] is None:
                minp[loc] = tallyval
                maxp[loc] = tallyval
            else:
                minp[loc] = min(minp[loc], tallyval)
                maxp[loc] = max(maxp[loc], tallyval)
            
            # test for matching target category depends on whether string or numeric
            if isstrvar:
                val = case[1].rstrip()
            else:
                val = case[1]
            if val == predval:   # hoping that predval is an integer if numeric
                targetcount[loc] += wt
    finally:
        curs.CClose()
    
    # make lists of case count and target count in ntile order
    caseslist = sorted(cases.items())
    targetcountlist = [item[1] for item in sorted(targetcount.items())]
    sumtargetcount = sum(targetcountlist)
    minp = sorted(minp.items())
    maxp = sorted(maxp.items())
    
    # calculate ntile target proportions and cumulative probabilities
    cumkt, cumtargetcount, targetpropor, minprop, maxprop, cumgain, lift = docalc(ntilelist, caseslist, targetcountlist, sumtargetcount, minp, maxp)
    
    # produce the final table, keeping it as a dataset if charts are requested
    tableds = dotable(anycharts, predprob, depvar, targetlabel, ntilelist, minprop, cumtargetcount, maxprop, targetcountlist, caseslist, targetpropor, cumkt, lift, cumgain)
    if anycharts:
        spss.Submit(f"""OMSEND TAG="{tableds}".""")
        docharts(predprob, depvar, predval, tableds, mainds, gainchart, responsechart, liftchart, targetlabel)
        

def dotable(anycharts, predprob, depvar, targetlabel, ntilelist, minprop, cumtargetcount, maxprop, 
        targetcountlist, caseslist, targetpropor, cumkt, lift, cumgain):
    """ Calculate the ntile tables and return the OMS tablename"""
    
    tableds = "D" + str(random.random())
    if anycharts:
        spss.Submit(f"""DATASET DECLARE {tableds}.
OMS SELECT TABLES /IF SUBTYPES="ntileanalysis"
    /DESTINATION FORMAT=SAV OUTFILE={tableds}
    /TAG="{tableds}".""")
    from spss import CellText
    ###from spss import FormatSpec    
    spss.StartProcedure(_("Ntile Analysis"))
    pt = spss.BasePivotTable(_("Ntile Analysis"), "ntileanalysis")
    spss.AddProcedureFootnotes(_(f"Predicted probabilities from variable {predprob} for dependent variable {depvar}"))
    spss.AddProcedureFootnotes(_(f"Target Category: {targetlabel}"))
    rowdim = pt.Append(spss.Dimension.Place.row, _("Ntile"))
    coldim =  pt.Append(spss.Dimension.Place.column, _("Statistics"))
    cols = [_("Upper Bound"), _("Minimum Probability"), _("Maximum Probability"), _("Count"), _("Cumulative Count"),
            _("Target Responses"), _("Cumulative Target Responses"), 
            _("Target Response Rate %"), _("Gain %"), _("Lift %")]
    pt.SetCategories(coldim, [CellText.String(v) for v in cols])    
    pt.SetCategories(rowdim, [CellText.String(v) for v in range(len(ntilelist))])
    for i in range(len(ntilelist)):
        pt.SetCellsByRow(CellText.String(str(i+1)), [CellText.Number(v) for v in 
            [ntilelist[i], minprop[i], maxprop[i], caseslist[i][1], cumkt[i], targetcountlist[i], 
            cumtargetcount[i], targetpropor[i], cumgain[i], lift[i]]])
    spss.EndProcedure()
    return tableds

def docalc(ntilelist, caseslist, targetcountlist, sumtargetcount, minp, maxp):
    """ Calculate table statistics"""
    
    cumkt = []
    cumtargetcount = []
    targetpropor = []
    cumtargetpropor = []
    minprop = []
    maxprop = []
    cumgain = []
    lift = []
    for ntile in range(len(ntilelist)):
        if ntile == 0:
            cumkt.append(caseslist[ntile][1])
            cumtargetcount.append(targetcountlist[ntile])
        else:
            cumkt.append(caseslist[ntile][1] + cumkt[ntile -1])
            cumtargetcount.append(targetcountlist[ntile] + cumtargetcount[ntile - 1])
        targetpropor.append(targetcountlist[ntile] / max(caseslist[ntile][1], .1) * 100.)
        cumtargetpropor.append(cumtargetcount[ntile] / max(cumkt[ntile], .1) * 100.)
        cumgain.append(cumtargetcount[ntile] / max(sumtargetcount, .1) * 100.)
        lift.append(cumtargetcount[ntile] / max(cumkt[ntile], .1) * 100.)
        minprop.append(minp[ntile][1])
        maxprop.append(maxp[ntile][1])
    return cumkt, cumtargetcount, targetpropor, minprop, maxprop, cumgain, lift

def dontiles(ntiles, predprob):
    """Get the variable ntiles"""
    
    incr = 100./ntiles
    ntilelist = frange(incr, 100., incr)
    ntilestr = " ".join(str(item) for item in ntilelist)
    cmd = f"""FREQ VARIABLES={predprob} /FORMAT=NOTABLE 
/PERCENTILES={ntilestr}."""
    tag, errlevel = spssaux.createXmlOutput(cmd, omsid='Frequencies', subtype="Statistics") # can't fail
    ntileupperbnds = spss.EvaluateXPath(tag, "/", '//group[@text="Percentiles"]/category/cell[@number]/@number')
    ntileupperbnds = [float(item) for item in ntileupperbnds]
    spss.DeleteXPathHandle(tag)
    if len(ntileupperbnds) <= 1:
        raise ValueError(_(f"There are fewer than two nonmissing cases.  Procedure stopped."))
    ntileupperbnds[-1] = 1.    # since these are floats, make sure the last value is exactly 1.
    nbnds = len(ntileupperbnds)
    return ntilelist, ntileupperbnds, nbnds
        
    
    
def frange(start, stop, step):
    "range with fractional values"
    res = []
    while start <= stop + 1e-10:
        res.append(float(start))
        start += step
    return res
    
def docharts(predprob, depvar, predval, tableds, mainds, gainchart, responsechart, liftchart, targetlabel):
    """Produce requested charts and close charting dataset
    
    predprob and depvar are the tabled variables
    tableds is the name of the charting dataset
    mainds is the name of the main cumulative
    gainchart, responsechart, and liftschart are the chart types"""
    
    # table structure variables
    # Command_,    Subtype_,    Label_,    
    # Var1 - ntile number
    #UpperBound
    #MinimumProbability
    #MaximumProbability
    #Count
    #CumulativeCount
    #TargetResponses
    #CumulativeTargetResponses
    #TargetResponseRate
    #Gain
    #Lift    
    
    gainchartcmd = f"""GGRAPH
  /GRAPHDATASET NAME="graphdataset" VARIABLES=UpperBound MEAN(Count)[name="Count"] MEAN(Gain)[name="MEAN_Gain"] 
  /GRAPHSPEC SOURCE=INLINE.
BEGIN GPL
  PAGE: begin(scale(500px, 400px))
  SOURCE: s=userSource(id("graphdataset"))
  DATA: UpperBound=col(source(s), name("UpperBound"))
  DATA: Count = col(source(s), name("Count"))
  DATA: MEAN_Gain=col(source(s), name("MEAN_Gain"))
  DATA: x = iter(0, 100, .1)
  TRANS: y = eval(x)
  ELEMENT: line(position(x*y), size(size.tiny), color(color.lightgray))
  GUIDE: axis(dim(1), label("Ntile Upper Bounds"))
  GUIDE: axis(dim(2), label("Gain %"))
  GUIDE: text.title(label("Gain for {predprob} Dependent Variable: {depvar}"))
  GUIDE: text.subtitle(label("Category: {targetlabel}"))
  GUIDE: legend(aesthetic(aesthetic.color.brightness), null())
  GUIDE: text.footnote(label("Darker points indicate smaller counts"))
  SCALE: linear(dim(2), include(0))
  ELEMENT: line(position(UpperBound*MEAN_Gain), missing.wings())
  ELEMENT: point(position(UpperBound*MEAN_Gain), size(size.medium), color(color.white), color.brightness(Count))  
  PAGE: end()
  END GPL."""
    
    responsechartcmd = f"""GGRAPH
  /GRAPHDATASET NAME="graphdataset" VARIABLES=UpperBound Count MEAN(TargetResponseRate)[name="TargetResponseRate"] 
  /GRAPHSPEC SOURCE=INLINE.
BEGIN GPL
  PAGE: begin(scale(500px, 400px))
  SOURCE: s=userSource(id("graphdataset"))
  DATA: UpperBound=col(source(s), name("UpperBound"), unit.category())
  DATA: TargetResponseRate=col(source(s), name("TargetResponseRate"))
  DATA: Count=col(source(s), name("Count"))
  GUIDE: axis(dim(1), label("Ntile Upper Bounds"))
  GUIDE: axis(dim(2), label("Target Response Rate %"))
  GUIDE: legend(aesthetic(aesthetic.size), null())
  GUIDE: text.title(label(" Response Rate for {predprob} Dependent Variable: {depvar}"))
  GUIDE: text.subtitle(label("Category: {targetlabel}"))
  GUIDE: text.footnote(label("Narrower bars indicate smaller counts"))
  SCALE: linear(dim(2), include(0))
  ELEMENT: interval(position(UpperBound*TargetResponseRate), size(Count))
  PAGE: end()
  END GPL."""
    liftchartcmd = f"""GGRAPH
  /GRAPHDATASET NAME="graphdataset" VARIABLES=UpperBound Count MEAN(Lift)[name="MEAN_lift"] 
  /GRAPHSPEC SOURCE=INLINE.
BEGIN GPL
  PAGE: begin(scale(500px, 400px))
  SOURCE: s=userSource(id("graphdataset"))
  DATA: UpperBound=col(source(s), name("UpperBound"))
  DATA: MEAN_lift=col(source(s), name("MEAN_lift"))
  DATA: Count=col(source(s), name("Count"))
  GUIDE: axis(dim(1), label("Ntile Upper Bounds"))
  GUIDE: axis(dim(2), label("Lift %"))
  GUIDE: text.title(label("Lift for {predprob} Dependent Variable: {depvar}"))
  GUIDE: text.subtitle(label("Category: {targetlabel}"))
  GUIDE: text.footnote(label("Darker points indicate smaller counts"))
  GUIDE: legend(aesthetic(aesthetic.color.brightness), null())
  SCALE: linear(dim(2), include(0))
  ELEMENT: line(position(UpperBound*MEAN_lift), missing.wings())
  ELEMENT: point(position(UpperBound*MEAN_lift), size(size.medium), color(color.white), color.brightness(Count))
  PAGE: end()
  END GPL."""

    spss.Submit(f"DATASET ACTIVATE {tableds}.")
    try:
        if gainchart:
            spss.Submit(gainchartcmd)
        if responsechart:
            spss.Submit(responsechartcmd)
        if liftchart:
            spss.Submit(liftchartcmd)
    finally:
        spss.Submit(f"""DATASET CLOSE {tableds}""")
        spss.Submit(f"""DATASET ACTIVATE {mainds}""")

    
def  Run(args):
    """Execute the STATS NTILE ANALYSIS command"""

    args = args[list(args.keys())[0]]


    oobj = Syntax([
        Template("PREDPROB", subc="",  ktype="existingvarlist", var="predprob", islist=False),
        Template("DEPVAR", subc="", ktype="existingvarlist", var="depvar"),
        Template("NTILES", subc="", ktype="int", var="ntiles",
            vallist=(2,)),
        Template("PREDVAL", subc="", ktype="str", var="predval"),
    
        Template("GAIN", subc="CHARTS", ktype="bool", var="gainchart"),
        Template("RESPONSERATE", subc="CHARTS", ktype="bool", var="responsechart"),
        Template("LIFT", subc="CHARTS", ktype="bool", var="liftchart")
    ])
        
    #debugging
    try:
        import wingdbstub
        if wingdbstub.debugger != None:
            import time
            wingdbstub.debugger.StopDebug()
            time.sleep(2)
            wingdbstub.debugger.StartDebug()
    except:
        pass

    #enable localization
    global _
    try:
        _("---")
    except:
        def _(msg):
            return msg

    # A HELP subcommand overrides all else
    if "HELP" in args:
        #print helptext
        helper()
    else:
        processcmd(oobj, args, dontile, vardict=spssaux.VariableDict())

def helper():
    """open html help in default browser window
    
    The location is computed from the current module name"""
    
    import webbrowser, os.path
    
    path = os.path.splitext(__file__)[0]
    helpspec = "file://" + path + os.path.sep + \
         "markdown.html"
    
    # webbrowser.open seems not to work well
    browser = webbrowser.get()
    if not browser.open_new(helpspec):
        print(("Help file not found:" + helpspec))
try:    #override
    from extension import helper
except:
    pass        
