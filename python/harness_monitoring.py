#!/bin/python
import os
import glob
import math
from array import array
import sys
import time
import subprocess
from optparse import OptionParser
import time
import datetime
import findFiles
import harness_definition

#print date
print("----------------------------------------------------------------------------------")
print(datetime.datetime.now())
print("----------------------------------------------------------------------------------")

start = time.time()

#parameters
current_dir = os.getcwd();

ignored_ntuples_label_list = ["obsolete","failed"]#ntuples containing anywhere in the path these labels will be ignored (eg corrupted files within the given tag_list)

#parse arguments
parser = OptionParser()
parser.add_option('--submit',          action='store_true',             dest='submit',          default=False,      help='submit jobs')
parser.add_option("--ntuple",          action="store",      type="str", dest="ntuple_dir",                          help="ntuple directory to run manually on")  
parser.add_option("--nfiles",          action="store",      type="int", dest="nfiles",                              help="ntuples to run on ")  
parser.add_option("--dbname",          action="store",      type="str", dest="dbname",                              help="dbname:: to retrieve files from ecalutomation, see doc")  
parser.add_option("--campaign",        action="store",      type="str", dest="campaign",                            help="campaing:: to retrieve files from ecalutomation, see doc")  
parser.add_option("--eras",            action="store",      type="str", dest="eras",                                help="eras:: to retrieve files from ecalutomation, see doc")  
parser.add_option("-q", "--queue",     action="store",      type="str", dest="condor_queue",    default="workday",  help="condor queue: espresso, longlunch, workday...")  
parser.add_option("-l", "--label",     action="store",      type="str", dest="label",                               help="job label")
parser.add_option("-v", "--verbosity", action="store",      type="int", dest="verbosity",       default=1,          help="verbosity level")
parser.add_option("-o", "--outdir",    action="store",      type="str", dest="outdir",          default="./",       help="output directory")
parser.add_option("-e", "--exedir",    action="store",      type="str", dest="exedir",          default="./bin/",   help="executable directory")
parser.add_option("-c", "--cfg",       action="store",      type="str", dest="configFile",                          help="template config file")
parser.add_option("-t", "--task",      action="store",      type="str", dest="tasklist",        default="runDivide,scaleMonitor",  
                 help="tasks to do it accepts buildTemplate,runDivide,scaleMonitor default:runDivide,scaleMonitor")
parser.add_option('--EE',              action='store_true',             dest='EE',              default=False,      help='run endcap calibration')
parser.add_option('--tier0',           action='store_true',             dest='tier0',           default=False,      help='submit to CAF queues (only if you are logged in lxplus-t0.cern.ch)')
parser.add_option('--groupByEras',     action='store_true',             dest='groupByEras',     default=False,      help='group files by Eras, eg Run2018C, Run2018D...')
parser.add_option('--groupByN',        action='store_true',             dest='groupByN',        default=False,      
                 help='group files in batches of a certain number, specify the number with --nfiles')
parser.add_option('--runManually',     action='store_true',             dest='runManually',     default=False,      
                 help='to run ntuples specified by a path provided, remember to specify the path with option --ntuple')
(options, args) = parser.parse_args()


#create outdir
eras = []
os.system("mkdir -p "+str(options.outdir))
if options.eras != None: eras = options.eras.split(',')

#get ecalelf ntuples for the calibration
if options.runManually:
    if options.ntuple_dir == None: 
        print ("Please, provide the ntuples path with --ntuple --> EXIT")
        sys.exit()

    selected_filelist, extracalibtree_filelist = findFiles.findFiles(options.ntuple_dir, "unmerged", ignored_ntuples_label_list)

    print (">>>>> Running on ecalelf ntuples manually produced: ")
    print (options.ntuple_dir)
 

    #Group the files together, NB when building templates you need statistics (at least an era)
    if options.groupByEras:
        if options.eras == None: 
            print ("Please, provide the eras with --eras --> EXIT")
            sys.exit()
        tag_list = eras 
        selected_filelist, extracalibtree_filelist = findFiles.groupFilesByTag(selected_filelist, extracalibtree_filelist, tag_list )

    elif options.groupByN:
        if options.nfiles == None: 
            print ("Please, insert the number of files with --nfiles --> EXIT")
            sys.exit() 
        selected_filelist, extracalibtree_filelist = findFiles.groupFiles(selected_filelist, extracalibtree_filelist, options.nfiles )

    else:  selected_filelist, extracalibtree_filelist = findFiles.groupFiles(selected_filelist, extracalibtree_filelist, len(selected_filelist) )

   
else: #getting from ecalautomation 
    print (">>>>> Running on ecalautomation ecalelf ntuples: ")
    print ("dbname            :   campaign    :                  eras                   :   grouping ")
    print (options.dbname, options.campaign, eras,  options.groupByEras)
    #Find files and leave them separated by ERA if specified (groupByEras == True), NB when building templates you need statistics (at least an era)
    selected_filelist, extracalibtree_filelist = findFiles.findFilesAuto(options.dbname, options.campaign, eras, options.groupByEras)

if (len(selected_filelist)>0):
    print("")
    if(options.verbosity>=1):
        print("-----------------------")
        print("ntuples fileList")
        for filename in selected_filelist:
            print (filename) 
        print("-----------------------")
        print("extraCalibTree filelist")
        for filename in extracalibtree_filelist:
            print (filename) 
        print("-----------------------")

else:
    print("")
    print("NOT any file found --> EXIT")
    sys.exit()

if(options.verbosity>=1):
    print ("grouped files")
    for filename in selected_filelist:
        print(filename) 

#create folder for the job
job_parent_folder=current_dir+"/jobs/"+str(options.label)+"/"
os.system("mkdir -p "+job_parent_folder)

#create the log folder
os.system("mkdir -p "+job_parent_folder+"/log/")

#get the harness ranges (when building templates, you need one template for each module, hence ieta varying while iphi 1-360)
if 'buildTemplate' in options.tasklist: harness_ranges =  harness_definition.GetModuleRanges()
else: harness_ranges = harness_definition.GetHarnessRanges()

print(">>>>>> Generating jobs for task lists "+options.tasklist)

#make the monitoring files .cfg, .sh, and .sub
for iFile in range(0,len(selected_filelist)):
    selected_filename=selected_filelist[iFile]
    extracalibtree_filename=extracalibtree_filelist[iFile]
    if(options.verbosity>=1):
        print(">>> Generating job for file "+selected_filename)
    for harness_range in harness_ranges:
        etamin = harness_range[0]
        etamax = harness_range[1]
        phimin = harness_range[2]
        phimax = harness_range[3]

        if(options.verbosity>=1):
            print(">>> Generating job for harness IEta_%i_%i_IPhi_%i_%i"%(etamin,etamax,phimin,phimax))

        jobdir="%s/IEta_%i_%i_IPhi_%i_%i/job_file_%i/"%(job_parent_folder,etamin,etamax,phimin,phimax,iFile)
        os.system("mkdir -p "+jobdir)
        outdir = "%s/IEta_%i_%i_IPhi_%i_%i/"%(options.outdir,etamin,etamax,phimin,phimax) 
        os.system("mkdir -p "+outdir)

        with open(str(options.configFile)) as fi:
            contents = fi.read()
            replaced_contents = contents.replace("SELECTED_INPUTFILE", selected_filename).replace("EXTRACALIBTREE_INPUTFILE", extracalibtree_filename)
            replaced_contents = replaced_contents.replace("IETAMIN",str(etamin)) 
            replaced_contents = replaced_contents.replace("IETAMAX",str(etamax))
            replaced_contents = replaced_contents.replace("IPHIMIN",str(phimin)) 
            replaced_contents = replaced_contents.replace("IPHIMAX",str(phimax)) 
            replaced_contents = replaced_contents.replace("OUTPUT_RUNDIVIDE","%s/out_file_%i_runranges.root"%(outdir,iFile))
            replaced_contents = replaced_contents.replace("OUTPUT_SCALEMONITORING","%s/out_file_%i_scalemonitoring.root"%(outdir,iFile))
            replaced_contents = replaced_contents.replace("OUTPUT_FOLDER",outdir)
            cfgfilename=jobdir+"/config.cfg"
            with open(cfgfilename, "w") as fo:
                fo.write(replaced_contents)

        ##### creates script #######
        outScriptName="%s/job_file_%i.sh"%(jobdir, iFile)
        outScript = open(outScriptName,"w")
        outScript.write("#!/bin/bash\n")
        #outScript.write('source setup.sh\n') 
        outScript.write("cd /afs/cern.ch/work/f/fcetorel/private/work2/EFlow/CMSSW_10_5_0/src/\n") #temporary solution
        outScript.write('eval `scram runtime -sh`\n');
        outScript.write("cd -\n");
        outScript.write("echo $PWD\n");
        outScript.write("%s/LaserMonitoring.exe --cfg %s"%(str(options.exedir),cfgfilename))
        for task in options.tasklist.split(','):
            outScript.write(" --"+task)
        outScript.write("\n")
        outScript.write("echo finish\n") 
        outScript.close();
        os.system("chmod 777 "+outScriptName)

#generate condor multijob submitfile for each task
condorsubFilename=job_parent_folder+"/submit_jobs.sub"
condorsub = open( condorsubFilename,"w")
condorsub.write("executable            = $(scriptname)\n")
condorsub.write("output                = $(scriptname).$(ClusterId).out\n")
condorsub.write("error                 = $(scriptname).$(ClusterId).err\n")
condorsub.write("log                   = %s/log/log.$(ClusterId).log\n"%job_parent_folder)
condorsub.write('+JobFlavour           = "%s"\n'%options.condor_queue)
if options.tier0:
    condorsub.write('+AccountingGroup      = "group_u_CMS.CAF.ALCA"\n')
condorsub.write("queue scriptname matching %s/IEta_*_*_IPhi_*_*/job_file_*/*.sh\n"%job_parent_folder)
condorsub.close()

submit_command = "condor_submit "+condorsubFilename
print("SUBMIT COMMAND: "+submit_command)

#submit in case the option is given
if(options.submit):
    os.system(submit_command)


end = time.time()

#print ("How much does it takes???", (end - start) , "s")

