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


def findFiles(ntuple_dir,ntuples_type,tag_list,ignored_ntuples_label_list,selected_filelist=[],extracalibtree_filelist=[]):
    #get ntuples for the calibration
    for root, dirs, files in os.walk(ntuple_dir):
        for file in files:
            if type_match(file,ntuples_type):
                if(not any (ignored_ntuples_label in os.path.join(root, file) for ignored_ntuples_label in ignored_ntuples_label_list)):
                    if file.find("extraCalibTree")==-1:
                        selected_filelist.append(os.path.join(root, file))
                        extracalibtree_filename = generate_extracalibtree_filename(file,ntuples_type)
                        extracalibtree_filelist.append(os.path.join(root, extracalibtree_filename))
    return selected_filelist,extracalibtree_filelist


def type_match(filename,ntuples_type):
    if not filename.endswith(".root"):
        return False
    if ntuples_type=="EGamma":
        if filename.startswith("EGamma"):
            return True
        else:
            return False
    if ntuples_type=="before2018":
        if filename.startswith("DoubleEG") or filename.startswith("SingleElectron"):
            return True
        else:
            return False
    if ntuples_type=="unmerged":
        if filename.startswith("ntuple"):
            return True

        else:
            return False

    return False

def generate_extracalibtree_filename(filename,ntuples_type):
    filenamecopy = filename
    if ntuples_type=="unmerged":
        if filename.startswith("ntuple_"): return filenamecopy.replace("ntuple_","extraCalibTree_")
        if filename.startswith("ntuple-"): return filenamecopy.replace("ntuple-","extraCalibTree-")
    else:
        return ("extraCalibTree-"+filenamecopy)

def groupFiles(selected_filelist, extracalibtree_filelist, Nfiles_per_group):
    print "grouping "+str(len(selected_filelist))+" files in groups of "+str(Nfiles_per_group)+" units"
    grouped_selected_filelist = []
    grouped_extracalibtree_filelist = []
    
    Nfiles_in_group=0
    selected_filename_group_str=""
    extracalibtree_filename_group_str=""
    for ifile in range(0,len(selected_filelist)):

        if Nfiles_in_group<Nfiles_per_group-1:#append file
            selected_filename_group_str += selected_filelist[ifile]+" \\ \n"
            extracalibtree_filename_group_str += extracalibtree_filelist[ifile]+" \\ \n"
            Nfiles_in_group=Nfiles_in_group+1

        elif Nfiles_in_group==Nfiles_per_group-1:#append file, close, and reset
                selected_filename_group_str += selected_filelist[ifile]
                extracalibtree_filename_group_str += extracalibtree_filelist[ifile]
                grouped_selected_filelist.append(selected_filename_group_str)
                grouped_extracalibtree_filelist.append(extracalibtree_filename_group_str)
                Nfiles_in_group=0
                selected_filename_group_str=""
                extracalibtree_filename_group_str=""

    print str(len(grouped_selected_filelist))+" groups created"
    return grouped_selected_filelist, grouped_extracalibtree_filelist

def groupFilesByTag(selected_filelist, extracalibtree_filelist, tag_list):
    print "grouping "+str(len(selected_filelist))+" files in groups sharing the same value among the following "
    print tag_list
    grouped_selected_filelist = []
    grouped_extracalibtree_filelist = []
    
    for tag in tag_list:
        Nfiles_in_group=0
        selected_filename_group_str=""
        extracalibtree_filename_group_str=""
        for ifile in range(0,len(selected_filelist)):
            if tag in selected_filelist[ifile]:
                selected_filename_group_str += selected_filelist[ifile]+" \\ \n"
                extracalibtree_filename_group_str += extracalibtree_filelist[ifile]+" \\ \n"
                Nfiles_in_group=Nfiles_in_group+1
        if Nfiles_in_group>0: 
            selected_filename_group_str = selected_filename_group_str[:-4]
            extracalibtree_filename_group_str = extracalibtree_filename_group_str[:-4]
            grouped_selected_filelist.append(selected_filename_group_str)
            grouped_extracalibtree_filelist.append(extracalibtree_filename_group_str)
        else:
            print "[WARNING]: can't find any file matching "+tag+" tag"
    print str(len(grouped_selected_filelist))+" groups created"
    return grouped_selected_filelist, grouped_extracalibtree_filelist

    

