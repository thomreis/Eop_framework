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

#print date
print("----------------------------------------------------------------------------------")
print(datetime.datetime.now())
print("----------------------------------------------------------------------------------")

#parameters
current_dir = os.getcwd();
ntuple_dir = "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/ecalelf/ntuples/13TeV/ALCARERECO/103X_dataRun2_v6_ULBaseForICs_newRegV1/"#parent folder containing all the ntuples of interest
#ntuple_dir = "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/ecalelf/ntuples/13TeV/ALCARERECO/102X_dataRun2_Sep2018Rereco_harnessCorr_newReg/"
#ntuple_dir="/home/fabio/work/Eop_framework/data/"
tag_list = ["Run2016C","Run2016D","Run2016E","Run2016F","Run2016G","Run2016H"]#tag for the monitoring = any label in the ntuple path identifying univoquely the ntuples of interest
#tag_list = ["Run2017C"] #tag for the monitoring
ignored_ntuples_label_list = ["obsolete"]#ntuples containing anywhere in the path these labels will be ignored (eg ntuples within a tag containing some error)


#parse arguments
parser = OptionParser()
parser.add_option('--submit',          action='store_true',             dest='submit',          default=False,      help='submit jobs')
parser.add_option("-l", "--label",     action="store",      type="str", dest="label",                               help="job label")
parser.add_option("-v", "--verbosity", action="store",      type="int", dest="verbosity",       default=1,          help="verbosity level")
parser.add_option("-o", "--outdir",    action="store",      type="str", dest="outdir",          default="./",       help="output directory")
parser.add_option(      "--jobdir",    action="store",      type="str", dest="jobdir",          default="",         help="job directory, if not provided use current directory")
parser.add_option("-e", "--exedir",    action="store",      type="str", dest="exedir",          default="./build/", help="executable directory")
parser.add_option("-c", "--cfg",       action="store",      type="str", dest="configFile",                          help="template config file")
parser.add_option("-N", "--Nloop",     action="store",      type="int", dest="Nloop",           default=15,         help="number of loop")
parser.add_option("--RestartFromLoop", action="store",      type="int", dest="RestartFromLoop", default=0,          help="restart existing calibration from the given loop")
parser.add_option('--odd',             action='store_true',             dest='odd',             default=False,      help='run only on odd entries')
parser.add_option('--even',            action='store_true',             dest='even',            default=False,      help='run only on even entries')
parser.add_option('--EE',              action='store_true',             dest='EE',              default=False,      help='run endcap calibration')
parser.add_option('--tier0',           action='store_true',             dest='tier0',           default=False,      help='submit to CAF queues (only if you are logged in lxplus-t0.cern.ch)')

(options, args) = parser.parse_args()

splitstat = ["odd","even"]
if(options.odd):
    splitstat = ["odd"]
if(options.even):
    splitstat = ["even"]

tasklist = ["BuildEopEta","ComputeIC"]

additional_options = ""
if options.EE:
    additional_options += " --EE "
else:
    print("setting up barrel calibration, if you want endcap calibration add the option --EE")

#create outdir
os.system("mkdir -p "+str(options.outdir))

#get ntuples for the calibration
selected_filelist,extracalibtree_filelist = findFiles.findFiles(ntuple_dir,"unmerged",tag_list,ignored_ntuples_label_list)

if (len(selected_filelist)>0):
    print()
    print("Run calibration on "+str(len(selected_filelist))+" files:")
    if(options.verbosity>=1):
        print("-----------------------")
        for filename in selected_filelist:
            print(filename)
        print("-----------------------")
        print("auto-generated extraCalibTree filelist")
        for filename in extracalibtree_filelist:
            print(filename) 
        print("-----------------------")

else:
    print()
    print("NOT any file found --> EXIT")
    sys.exit()

#if running on unmerged files i need to reduce the number of jobs to submit to a reasonable value --> Group the files together
if len(selected_filelist)>200:
    selected_filelist, extracalibtree_filelist = findFiles.groupFiles(selected_filelist, extracalibtree_filelist, int(len(selected_filelist)/25) )
    if(options.verbosity>=1):
        print("grouped files")
        for filename in selected_filelist:
            print(filename) 

#create folder for the job
print(options)
if str(options.jobdir)=="":
    job_parent_folder=current_dir+"/jobs/"+str(options.label)+"/"
else:
    job_parent_folder=options.jobdir+"/jobs/"+str(options.label)+"/"
os.system("mkdir -p "+job_parent_folder)

#create the log folder
os.system("mkdir -p "+job_parent_folder+"/log/")

#create DAGMan file to manage submitting hierarchy
dagFilename=job_parent_folder+"/submit_manager.dag"
dagFile = open( dagFilename,"w")

#make the monitoring files .cfg, .sh, and .sub
for iLoop in range(options.RestartFromLoop,options.Nloop):
    print("> Generating job for loop "+str(iLoop))

    for task in tasklist:
        if(options.verbosity>=1):
            print(">> Generating job for "+task)

        for iFile in range(0,len(selected_filelist)):
            selected_filename=selected_filelist[iFile]
            extracalibtree_filename=extracalibtree_filelist[iFile]
            if(options.verbosity>=1):
                print(">>> Generating job for file "+selected_filename)
                
            jobdir=job_parent_folder+"/job_loop_"+str(iLoop)+"_file_"+str(iFile)+"_"+task+"/"
            os.system("mkdir "+jobdir)
            
            with open(str(options.configFile)) as fi:
                contents = fi.read()
                replaced_contents = contents.replace("SELECTED_INPUTFILE", selected_filename).replace("EXTRACALIBTREE_INPUTFILE", extracalibtree_filename)
            cfgfilename=jobdir+"/config.cfg"
            with open(cfgfilename, "w") as fo:
                fo.write(replaced_contents)

            for split in splitstat:
                ##### creates executable options #######
                BUILDEOPETA_OUTPUT= str(options.outdir)+"/EopEta_loop_"+str(iLoop)+"_file_"+str(iFile)+"_"+split+".root"
                UPDATEIC_OUTPUT= str(options.outdir)+"/IC_loop_"+str(iLoop)+"_file_"+str(iFile)+"_"+split+".root"
                BUILDEOPETA_INPUT_OPTION=""
                UPDATEIC_INPUT_OPTION=""
                EOPWEIGHTRANGE_OPTION=""
                if iLoop==0:
                    if "BuildEopEta" in task:
                        BUILDEOPETA_INPUT_OPTION=""
                        UPDATEIC_INPUT_OPTION=""
                        if not options.EE:
                            EOPWEIGHTRANGE_OPTION="--Eopweightrange 0.2 1.9 --Eopweightbins 204"
                        else:
                            EOPWEIGHTRANGE_OPTION="--Eopweightrange 0.1 3.0 --Eopweightbins 250"                        
                    if "ComputeIC" in task:
                        BUILDEOPETA_INPUT_OPTION="--Eopweight TH2F EopEta "+str(options.outdir)+"/EopEta_loop_"+str(iLoop)+".root"
                        UPDATEIC_INPUT_OPTION=""
                        EOPWEIGHTRANGE_OPTION=""
                else:
                    if "BuildEopEta" in task:
                        BUILDEOPETA_INPUT_OPTION="--Eopweight TH2F EopEta "+str(options.outdir)+"/EopEta_loop_"+str(iLoop-1)+".root"
                        UPDATEIC_INPUT_OPTION="--inputIC IC "+str(options.outdir)+"/IC_loop_"+str(iLoop-1)+".root"                    
                        if not options.EE:
                            EOPWEIGHTRANGE_OPTION="--Eopweightrange 0.85 1.15 --Eopweightbins 36"
                        else:
                            EOPWEIGHTRANGE_OPTION="--Eopweightrange 0.1 3.0 --Eopweightbins 250"                        

                    if "ComputeIC" in task:
                        BUILDEOPETA_INPUT_OPTION="--Eopweight TH2F EopEta "+str(options.outdir)+"/EopEta_loop_"+str(iLoop)+".root"
                        UPDATEIC_INPUT_OPTION="--inputIC IC "+str(options.outdir)+"/IC_loop_"+str(iLoop-1)+".root"                    
                        EOPWEIGHTRANGE_OPTION=""

                ##### creates script #######
                outScriptName=jobdir+"/job_file_"+str(iFile)+"_"+split+".sh"
                outScript = open(outScriptName,"w")
                outScript.write("#!/bin/bash\n")
                #outScript.write('source setup.sh\n')
                outScript.write("#cd /afs/cern.ch/work/f/fmonti/flashggNew/CMSSW_10_5_0/\n")
                outScript.write('#eval `scram runtime -sh`\n');
                outScript.write("#cd -\n");
                outScript.write("echo $PWD\n");
                outScript.write(
                    str(options.exedir)+"/"+task+".exe"+
                    " --cfg "+cfgfilename+
                    " "+UPDATEIC_INPUT_OPTION+
                    " "+BUILDEOPETA_INPUT_OPTION+
                    " --BuildEopEta_output "+BUILDEOPETA_OUTPUT+
                    " "+EOPWEIGHTRANGE_OPTION+
                    " --ComputeIC_output "+UPDATEIC_OUTPUT+
                    " --"+split+
                    " "+additional_options+"\n")
                outScript.write("echo finish\n") 
                outScript.close();
                os.system("chmod 777 "+outScriptName)

        #generate condor multijob submitfile for each task
        condorsubFilename=job_parent_folder+"/submit_"+task+"_loop_"+str(iLoop)+".sub"
        condorsub = open( condorsubFilename,"w")
        condorsub.write("executable            = $(scriptname)\n")
        condorsub.write("output                = $(scriptname).$(ClusterId).out\n")
        condorsub.write("error                 = $(scriptname).$(ClusterId).err\n")
        condorsub.write("log                   = "+job_parent_folder+"/log/log.$(ClusterId).log\n")
        condorsub.write('+JobFlavour           = "workday"\n')
        if options.tier0:
            condorsub.write('+AccountingGroup      = "group_u_CMS.CAF.ALCA"\n')
        condorsub.write("queue scriptname matching "+job_parent_folder+"/job_loop_"+str(iLoop)+"_file_*_"+task+"/*.sh\n")
        condorsub.close()
        #fill the submitting manager file
        dagFile.write("JOB "+task+"_loop_"+str(iLoop)+" "+condorsubFilename+"\n")

        #submit the merging step
        if "BuildEopEta" in task:
            mergescriptName=job_parent_folder+"/merge_"+task+"_loop_"+str(iLoop)+".sh"
            mergescript = open( mergescriptName,"w")
            mergescript.write("#!/bin/bash\n")
            mergescript.write("#cd /afs/cern.ch/work/f/fmonti/flashggNew/CMSSW_10_5_0/\n")
            mergescript.write('#eval `scram runtime -sh`\n');
            mergescript.write("#cd -\n");
            mergescript.write("hadd -f -k "+str(options.outdir)+"/EopEta_loop_"+str(iLoop)+".root "+str(options.outdir)+"/EopEta_loop_"+str(iLoop)+"_file_*_*.root\n")
            mergescript.write(str(options.exedir)+"/NormalizeBuildEopEta.exe --Eopweight TH2F EopEta "+str(options.outdir)+"/EopEta_loop_"+str(iLoop)+".root\n")
            mergescript.close()
            os.system("chmod 777 "+mergescriptName)
        if "ComputeIC" in task:
            mergescriptName=job_parent_folder+"/merge_"+task+"_loop_"+str(iLoop)+".sh"
            mergescript = open( mergescriptName,"w")
            mergescript.write("#!/bin/bash\n")
            mergescript.write("#cd /afs/cern.ch/work/f/fmonti/flashggNew/CMSSW_10_5_0/\n")
            mergescript.write('#eval `scram runtime -sh`\n');
            mergescript.write("#cd -\n");
            mergescript.write("hadd -f -k "+str(options.outdir)+"/IC_loop_"+str(iLoop)+".root "+str(options.outdir)+"/IC_loop_"+str(iLoop)+"_file_*_*.root\n")
            if iLoop==0:
                mergescript.write(str(options.exedir)+"/UpdateIC.exe --newIC IC "+str(options.outdir)+"/IC_loop_"+str(iLoop)+".root\n")
            else:
                mergescript.write(str(options.exedir)+"/UpdateIC.exe --oldIC IC "+str(options.outdir)+"/IC_loop_"+str(iLoop-1)+".root --newIC IC "+str(options.outdir)+"/IC_loop_"+str(iLoop)+".root\n")
            mergescript.close()
            os.system("chmod 777 "+mergescriptName)

        mergesubFilename=job_parent_folder+"/submit_merge"+task+"_loop_"+str(iLoop)+".sub"
        mergesub = open( mergesubFilename,"w")
        mergesub.write("executable            = "+mergescriptName+"\n")
        mergesub.write("output                = "+mergescriptName+".$(ClusterId).out\n")
        mergesub.write("error                 = "+mergescriptName+".$(ClusterId).err\n")
        mergesub.write("log                   = "+job_parent_folder+"/log/log.$(ClusterId).log\n")
        mergesub.write('+JobFlavour           = "longlunch"\n')
        if options.tier0:
            mergesub.write('+AccountingGroup      = "group_u_CMS.CAF.ALCA"\n')
        mergesub.write("queue 1\n")
        mergesub.close()

        #fill the submitting manager file
        dagFile.write("JOB merge"+task+"_loop_"+str(iLoop)+" "+mergesubFilename+"\n")
            
#setting hierarchy of the submitting manager file
for iLoop in range(options.RestartFromLoop,options.Nloop):
    for iTask in range(0,len(tasklist)):
        dagFile.write("PARENT "+tasklist[iTask]+"_loop_"+str(iLoop)+" CHILD merge"+tasklist[iTask]+"_loop_"+str(iLoop)+"\n")
        if iTask<(len(tasklist)-1):
            dagFile.write("PARENT merge"+tasklist[iTask]+"_loop_"+str(iLoop)+" CHILD "+tasklist[iTask+1]+"_loop_"+str(iLoop)+"\n")
        else:
            if( iLoop < (options.Nloop-1) ):
                dagFile.write("PARENT merge"+tasklist[iTask]+"_loop_"+str(iLoop)+" CHILD "+tasklist[0]+"_loop_"+str(iLoop+1)+"\n")

#add possibility to re-submit failed jobs
for iLoop in range(options.RestartFromLoop,options.Nloop):
    for iTask in range(0,len(tasklist)):
        dagFile.write("Retry "+tasklist[iTask]+"_loop_"+str(iLoop)+" 3\n")
        dagFile.write("Retry merge"+tasklist[iTask]+"_loop_"+str(iLoop)+" 3\n")

dagFile.close()

submit_command = "condor_submit_dag "+dagFilename
print("SUBMIT COMMAND: "+submit_command)
#submit in case the option is given
if(options.submit):
    os.system(submit_command)
