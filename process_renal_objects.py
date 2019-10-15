#!/usr/bin/python
##
# \file         process_renal_objects.py
# \author       Bill Hill
# \date         October 2019
# \version      $Id$
# \par
# Address:
#               MRC Human Genetics Unit,
#               MRC Institute of Genetics and Molecular Medicine,
#               University of Edinburgh,
#               Western General Hospital,
#               Edinburgh, EH4 2XU, UK.
# \par
# Copyright (C), [2019],
# The University Court of the University of Edinburgh,
# Old College, Edinburgh, UK.
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be
# useful but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free
# Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.
# \brief        Pre-processes, registers and combines renal vesicle
#               and s-shaped body assay images.
#               Woolz spatial domain object.
#               This script makes use of ANTs, many Woolz binaries
#               and PyWoolz, see:
#                 https://github.com/ANTsX/ANTs
#                 https://github.com/ma-tech/Woolz
#                 https://github.com/ma-tech/PyWoolz
##

from __future__ import print_function
import os
import re
import sys
import csv
import shutil
import commands
import traceback
import ctypes as c
import Wlz as w

libc = c.CDLL('libc.so.6')

class WlzError(Exception):
  pass

# Set these processing requirement for a run of the pipeline. In general
# those that come lower in the list require those above to have been run.
verbose               = True
dummy_run             = False
exit_on_cmd_error     = True
prepare_data          = False
compute_registration  = False
register_models       = False
check_thresholds_set  = False
create_gx_band_images = False
combine_to_model      = False
create_point_clouds   = False
create_ref_surface    = False
create_compound       = False
cut_sections          = False
compute_displacements = False
compute_gx_intersect  = False

assay_idx_range       = [0,131] # Assay range selection
model_idx_set         = [0,3]   # Model range selection

# Set base_dir to directory which contains the following directory
# hierarchy:
#   Registration/
#     Models/
#     Tmp/
#     Visualization/
#     Working/
#     ManualReg/
#     Thresholds/
# Note: You must set your base directory before running this script
base_dir              = /some/path

reg_dir               = base_dir + '/Registration'
mod_dir               = reg_dir + '/Models'
tmp_dir               = reg_dir + '/Tmp'
vis_dir               = reg_dir + '/Visualization'
work_dir              = reg_dir + '/Working'
lmk_dir               = reg_dir + '/ManualReg'
thresh_dir            = reg_dir + '/Thresholds'

# List of CSV files with manually defined threshold levels
thresh_files          = ['Threshold_RV_SSB_v2.csv', \
                         'Threshold_RV_SSB_v4.csv', \
                         'Threshold_RV_SSB_v2_01.csv' ]

# Temporary file base name
tmp_file              = tmp_dir + '/tmp'

# Number of threshold bands
n_thr_band            = 5

# Threshold band dilation
thr_band_dilation     = 2

# Cubic voxel size after resamppling
cubic_sz              = 0.34605

# Woolz affine transform required to achieve cubic voxels
# prompt% WlzFacts Tmp/cubic_tr.wlz
#   Object type: WLZ_AFFINE_TRANS.
#      Linkcount: 0.
#      Transform type: WLZ_TRANSFORM_3D_AFFINE.
#      Linkcount: 1.
#      mat: 0.547302   0          0          0         
#           0          0.547302   0          0         
#           0          0          0.988889   0         
#           0          0          0          1         
#      Property list NULL.
cubic_tr_file         = tmp_dir + '/cubic_tr.wlz'

# Template filename for manual registration landmark files as saved
# by WlzWarp
landmark_tmpl         = 'WlzWarp_landmarks_Image{:d}_to_Image{:d}.lmk'

# Enumeration of object types
class ObjType: #{
  unknown = 0
  rv      = 1         # Renal Vesicle
  ssb     = 2         # S-Shaped Body
#}

obj_type              = ['unknown', 'rv', 'ssb']

# Enumeration of image channels
class Channel: #{
  unknown = 0
  cdh1    = 1
  dapi    = 2
  jag1    = 3
  lef1    = 4
  sall1   = 5
  six2    = 6
  sox9    = 7
  troma1  = 8
  wt1     = 9
  lhx1    = 10
  pax2    = 11
  hnf1b   = 12
  foxc2   = 13
  cdh6    = 14
  cldn5   = 15
  emx2    = 16
  erbb4   = 17
  mafb    = 18
  mecom   = 19
  pappa2  = 20
  pou3f3  = 21
#}

# Reverse channel lookup
channel               = ['unknown' , 'cdh1'   , 'dapi'  , 'jag1'   , 'lef1'  , \
                         'sall1'   , 'six2'   , 'sox9'  , 'troma1' , 'wt1'   , \
                         'lhx1'    , 'pax2'   , 'hnf1b' , 'foxc2'  , 'cdh6'  , \
                         'cldn5'   , 'emx2'   , 'erbb4' , 'mafb'   , 'mecom' , \
                         'pappa2'  , 'pou3f3']

def matchString(s,a): #{
  match = None
  for i in range(0,len(a)): #{
    if s == a[i]: #{
      match = i
      break
    #}
  #}
  return match
#}

def matchObjType(s): #{
  return matchString(s, obj_type)
#}

def vrbMsg(s): #{
  if verbose: #{
    print(s)
  #}
#}

def landmarkFile(src, tgt): #{
  f = lmk_dir + '/' + landmark_tmpl.format(src, tgt)
  vrbMsg('checking for landmark file ' + f)
  if not os.path.isfile(f): #{
    f = None
  #}
  vrbMsg('landmark file exists = ' + str(bool(f)))
  return f
#}

def matchChannel(s): #{
  return matchString(s, channel)
#}

# Enumeration of species
class Species: #{
  unknown = 0
  human   = 1
  mouse   = 2
#}

species               = ['unknown', 'human', 'mouse']

def matchSpecies(s): #{
  return matchString(s, species)
#}

def doCmd(s): #{
  rtn = ''
  vrbMsg('doCmd s = ' + s)
  if not dummy_run: #{
    so = commands.getstatusoutput(cmd)
    vrbMsg('doCmd so = ' + str(so))
    if(not so[0] == 0): #{
      print('Error: ' + cmd, file=sys.stderr)
      if(exit_on_cmd_error): #{
        sys.exit(so[0])
      #}
    #}
    rtn = so[1]
  #}
  vrbMsg('doCmd rtn = ' + rtn)
  return(rtn)
#}

def copyFile(src, tgt): #{
  vrbMsg('copyFile ' + src + ' ' + tgt)
  if not dummy_run: #{
    shutil.copyfile(src, tgt)
  #}
#}

def renameFile(src, tgt): #{
  vrbMsg('renameFile ' + src + ' ' + tgt)
  if not dummy_run: #{
    os.rename(src, tgt)
  #}
#}

def orgFileNii(asy, img): #{
  nii_file = base_dir + '/' + asy['dir'] + '/' + img['orig_file']
  return nii_file
#}

def orgFileWlz(asy, img): #{
  nii_file = orgFileNii(asy, img)
  wlz_file = re.sub('\.nii$', '.wlz', nii_file)
  return wlz_file
#}

def drvFileBse(asy, img, drv): #{
  bse_file = work_dir + '/asy-' + ('%08d' % asy['index']) + '-' + drv + \
             '-' + channel[img['channel']]
  return bse_file
#}


def blrFileBse(asy, img): #{
  bse_file = drvFileBse(asy, img, 'blr')
  return bse_file
#}

def blrFileNii(asy, img): #{
  nii_file = blrFileBse(asy, img) + '.nii'
  return nii_file
#}

def blrFileWlz(asy, img): #{
  wlz_file = blrFileBse(asy, img) + '.wlz'
  return wlz_file
#}

def sclFileBse(asy, img): #{
  bse_file = drvFileBse(asy, img, 'scl')
  return bse_file
#}

def sclFileNii(asy, img): #{
  nii_file = sclFileBse(asy, img) + '.nii'
  return nii_file
#}

def sclFileWlz(asy, img): #{
  wlz_file = sclFileBse(asy, img) + '.wlz'
  return wlz_file
#}

def sclFileChnBase(asy, chn): #{
  scl_file = work_dir + '/asy-' + ('%08d' % asy['index']) + '-scl-' + \
             str(channel[chn])
  return scl_file
#}

def sclFileChnNii(asy, chn): #{
  scl_file = sclFileChnBase(asy, chn) + '.nii'
  return scl_file

def sclFileChnWlz(asy, chn): #{
  scl_file = sclFileChnBase(asy, chn) + '.wlz'
  return scl_file

def blrFileChnBase(asy, chn): #{
  blr_file = work_dir + '/asy-' + ('%08d' % asy['index']) + '-blr-' + \
             str(channel[chn])
  return blr_file
#}

def blrFileChnWlz(asy, chn, must_exist): #{
  wlz_file = blrFileChnBase(asy, chn) + '.wlz'
  if bool(must_exist) and (not os.path.exists(wlz_file)): #{
    wlz_file = None
  #}
  return wlz_file
#}

def blrFileChnNii(asy, chn): #{
  nii_file = blrFileChnBase(asy, chn) + '.nii'
  return nii_file
#}

def dspFile(mod, asy): #{
  dsp_file = work_dir + '/dsp-' + str(mod['index']) + '-' + \
             str(mod['assays'][0]) + '-' + str(asy['index']) + \
             '.wlz'
  return dsp_file
#}

def wlzAfnTfmFile(mod, asy): #{
  wlz_file = work_dir + '/afn-' + str(mod['index']) + '-' + \
             str(mod['assays'][0]) + '-' + str(asy['index']) + '.wlz'
  return wlz_file
#}

def tfmFile(mod, asy): #{
  nii_file = work_dir + '/tfm-' + str(mod['index']) + '-' + \
             str(mod['assays'][0]) + '-' + str(asy['index']) + \
             '.nii'
  return nii_file
#}

def regFileChnBase(asy, chn, blr): #{
  reg_file = work_dir + '/reg-' + ('%08d' % asy['index'])
  if bool(blr): #{
    reg_file = reg_file + '-blr-'
  else: #}{
    reg_file = reg_file + '-scl-'
  #}
  reg_file = reg_file + str(channel[chn])
  return reg_file
#}

def regFileChnNii(asy, chn, blr): #{
  nii_file = regFileChnBase(asy, chn, blr) + '.nii'
  return nii_file
#}

def regFileChnWlz(asy, chn, blr): #{
  wlz_file = regFileChnBase(asy, chn, blr) + '.wlz'
  return wlz_file
#}

def secFileChnTif(asy_idx, angles, chn): #{
  tif_file = work_dir + '/sec-' + ('%03d' % asy_idx) + '-' + angles + '-' + \
             str(channel[chn]) + '.tif'
  return tif_file
#}

def regPtsFileChnVtk(asy, chn_idx): #{
  vtk_file = work_dir + '/reg-pts-' + ('%08d' % asy['index']) + '-' + \
             channel[chn_idx] + '.vtk'
  return(vtk_file)
#}

def modFileBlrChnWlz(mod_idx, chn_idx): #{
  wlz_file = work_dir + '/mod-blr-' + str(mod_idx) + '-' + channel[chn_idx] + \
             '.wlz'
  return wlz_file
#}

def modFileHeqChnWlz(mod_idx, chn_idx): #{
  wlz_file = work_dir + '/mod-heq-' + str(mod_idx) + '-' + channel[chn_idx] + \
             '.wlz'
  return wlz_file
#}

def modFileThrChnWlz(mod_idx, chn_idx, blr, typ): #{
  typ = int(typ)
  assert typ in range(0,3)
  nam = 'thr'
  val = ''
  if bool(blr): #{
    nam = 'thrblr'
  #}
  if typ == 0: #{
    pass
  elif typ == 1: #}{
    nam = nam + 'max'
  elif typ == 2: #}{
    nam = 'occ'
    val = '-1'
  #}
  wlz_file = work_dir + '/mod-' + nam + '-' + str(mod_idx) + '-' + \
             channel[chn_idx] + val + '.wlz'
  return wlz_file
#}

def modFileSclChnWlz(mod_idx, chn_idx): #{
  wlz_file = work_dir + '/mod-scl-' + str(mod_idx) + '-' + channel[chn_idx] + \
             '.wlz'
  return wlz_file
#}

def modPtsFileChnVtk(mod_idx, chn_idx): #{
  vtk_file = work_dir + '/mod-pts-' + str(mod_idx) + '-' + channel[chn_idx] + \
             '.vtk'
  return(vtk_file)
#}

def modSrfFileChnVtk(mod_idx, chn_idx): #{
  vtk_file = work_dir + '/mod-srf-' + str(mod_idx) + '-' + channel[chn_idx] + \
             '.vtk'
  return(vtk_file)
#}

def modCpdFile(mod_idx): #{
  wlz_file = work_dir + '/mod-cpd-' + str(mod_idx) + '.wlz'
  return(wlz_file)
#}

def thrFileAsyChnWlz(asy_idx, chn_idx): #{
  wlz_file = work_dir + '/thr-' + ('%08d' % asy_idx) + channel[chn_idx] + \
             '.wlz'
  return(wlz_file)
#}

def modDspFileBody(mod_idx): #{
  wlz_file = work_dir + '/mod-dsp-' + str(mod_idx) + '-'
  return(wlz_file)
#}

def modFileThrOccWlz(mod_idx, chn_idx, thr_v): #{
  wlz_file = work_dir + '/mod-occ-' + str(mod_idx) + '-' + \
             channel[chn_idx] + '-' + str(thr_v) + '.wlz'
  return(wlz_file)
#}

def isnFile(mod_idx, typ): #{
  typ = int(typ)
  assert typ in range(0,3)
  t_str = ['', 'max_', 'occ_']
  isn_file = work_dir + '/intersections_' + t_str[typ] + str(mod_idx) + '.csv'
  return(isn_file)
#}

def assayByIndex(g_idx): #{
  idx = g_idx
  n = len(assays)
  if g_idx < 0: #{
    idx = 0
  elif g_idx >= n: #}{
    idx = n - 1
  #}
  asy = assays[idx]
  while (idx >= 0) and (asy['index'] > g_idx): #{
    idx = idx - 1
    asy = assays[idx]
  #}
  while (idx < n) and (asy['index'] < g_idx): #{
    idx = idx + 1
    asy = assays[idx]
  #}
  if not (asy['index'] == g_idx): #{
    asy = None
    print('Error: There is no assay with index ' + str(g_idx), file=sys.stderr)
  #}
  return(asy)
#}

# Get maximum value within the given Woolz object
def objMaxValue(obj): #{
  max = 0
  try: #{
    pv = [w.WlzPixelV() for i in range(2)]
    errNum = w.WlzGreyRange(obj, c.byref(pv[0]), c.byref(pv[1]))
    if bool(errNum): #{
      raise WlzError()
    #}
    # Here we know we'll only have integral value types
    max = {
      int(w.WLZ_GREY_INT):   lambda v: int(v.inv),
      int(w.WLZ_GREY_SHORT): lambda v: int(v.shv),
      int(w.WLZ_GREY_UBYTE): lambda v: int(v.ubv),
    }[pv[1].type](pv[1].v)
  #}
  except: #}{
    vrbMsg('Failed to find maximum object value.')
    if(exit_on_cmd_error): #{
      sys.exit(so[0])
    #}
  #}
  return(max)
#}

# Read assay data definitions and convert strings to 'enum' class members
execfile('assays.py')
for asy in assays: #{
  asy['species'] = matchSpecies(asy['species'])
  asy['obj_type'] = matchObjType(asy['obj_type'])
  for img in asy['images']: #{
    img['channel'] = matchChannel(img['channel'])
  #}
#}

# Read model definitions and convert strings to 'enum' class members
execfile('models.py')
for mod in models: #{
  mod['obj_type'] = matchObjType(mod['obj_type'])
  mod['species'] = matchSpecies(mod['species'])
#}

# Import thresholds from CSV file and put into assays
for tf in thresh_files: #{
  thresh_file = thresh_dir + '/' + tf
  f = open(thresh_file, 'rt')
  try: #{
    line = 1
    thresh =  csv.reader(f)
    for rec in thresh: #{
      # There are two possible threshold file formats here:
      #   v2: assay_idx, channel, value
      # or
      #   v4: assay_idx, original file, value, species, obj_type, 16 bit value
      if rec[0].isdigit(): #{
        match = False
        i = int(rec[0])
        asy = assayByIndex(i)
        g = None
        t = 0.0
        ts = '0.0'
        if (len(rec) == 3): #{
          # v2 threshold file
          g = ''.join(rec[1].lower().split())
          ts = rec[2].lower()
        elif (len(rec) == 6): #}{
          # v4 threshold file
          g = re.sub('^.*_([A-Z0-9]+)\....*$', '\\1', rec[1]).lower()
          if not (g == 'dapi'): #{
            ts = rec[2].lower()
          #}
        #}
        if ts == 'na' or ts == 'n/a': #{
          if not (g == 'dapi'): #{
            t = 255.0
          #}
        else: #}{
          t = float(ts)
        #}
        if bool(g): #{
          # For some reason we've used antibody troma1 not gene krt8
          if g == 'krt8': #{
            g = 'troma1'
          #}
          img = None
          for img in asy['images']: #{
            if channel[img['channel']] == g: #{
              img['thresh_8bit'] = t
              match = True
            #}
          #}
        #}
        vrbMsg('Threshold ' + str(asy['index'])  + ' ' + str(g) + ' ' + str(t))
        if not match: #{
          print('Error: No match in threshold file ' + thresh_file + \
                ', line ' + str(line) + ' ' + str(rec))
          sys.exit(1)
        #}
      #}
      line = line + 1
    #}
  finally: #}{
    f.close()
  #}
#}

# Check all thresholds are set
if check_thresholds_set: #{
  for asy in assays: #{
    vrbMsg('Assay index ' + str(asy['index']))
    for img in asy['images']: #{
      if not img['channel'] == Channel.dapi: #{
        if img['thresh_8bit'] == 0: #{
          print('Error: Threshold not set for assay ' + str(asy['index']) + \
                ', channel ' + channel[(img['channel'])])
          sys.exit(1)
        #}
      #}
    #}
  #}
#}

if prepare_data: #{
  vrbMsg('Begin Section - prepare_data')
  # Create output directories if they don't already exist
  cmd = 'mkdir -p ' + mod_dir + ' ' + tmp_dir + ' ' + vis_dir + ' ' + work_dir
  doCmd(cmd)
  # Convert the original files to Woolz format.
  # Sample the objects to get cubic voxels
  # Create sampled blurred DAPI and JAG1 images ready for registration.
  for asy in assays: #{
    vrbMsg('Assay index ' + str(asy['index']))
    voxRescale = None
    if asy['index'] >= assay_idx_range[0] and \
       asy['index'] <= assay_idx_range[1]: #{
      for img in asy['images']: #{
        vrbMsg('Original file ' + img['orig_file'])
        img_file_nii = orgFileNii(asy, img)
        img_file_wlz = orgFileWlz(asy, img)
        # Convert the NIfTI data to Woolz.
        cmd = 'WlzExtFFConvert  -f nii -F wlz ' + \
                              ' -o ' + img_file_wlz + ' ' + img_file_nii
        doCmd(cmd)
        if not bool(voxRescale): #{
          # Compute scale transform to give uniform voxels of required size,
          # with the size set by cubic_sz. Assume that we only need to do
          # this for the first image of the assay
          cmd = 'WlzFacts ' + img_file_wlz + ' 2>&1 | grep VoxelSize'
          voxSize = doCmd(cmd).split()[1:]
          vrbMsg('voxSize = ' + str(voxSize))
          voxRescale = [None] * len(voxSize)
          for i in range(0,len(voxSize)): #{
            voxRescale[i] = float(voxSize[i]) / cubic_sz
          #}
          vrbMsg('voxRescale = ' + str(voxRescale))
          tr_txt_file = tmp_file + '-pre-0.txt'
          cmd = 'printf \'' + str(voxRescale[0]) + ' 0 0 0\\n' + \
                              '0 ' + str(voxRescale[1]) + ' 0 0\\n' + \
                              '0 0 ' + str(voxRescale[2]) + ' 0\\n' + \
                              '0 0 0 1.0\\n\' > ' + tr_txt_file
          doCmd(cmd)
          cmd = ' WlzAffineTransformObj -3 -N -m ' + tr_txt_file + \
                ' -T ' + cubic_tr_file
          doCmd(cmd)
          os.remove(tr_txt_file)
        #}
        # Scale the data so that we have cubic voxels.
        scl_file_nii = sclFileNii(asy, img)
        scl_file_wlz = sclFileWlz(asy, img)
        cmd = 'WlzAffineTransformObj -t ' + cubic_tr_file + \
                  ' -o ' + scl_file_wlz + ' ' + img_file_wlz
        doCmd(cmd)
        # Convert the scaled image back to NIfTI
        cmd = 'WlzExtFFConvert -f wlz -F nii ' + \
                  ' -o ' + scl_file_nii + ' ' + scl_file_wlz
        doCmd(cmd)
        # Create blurred normalised images. Needed for registration and useful
        # for visualization.
        chn = img['channel']
        blr_file_wlz = blrFileWlz(asy, img)
        blr_file_nii = blrFileNii(asy, img)
        cmd = 'WlzSepFilterObj -m 2,2,2 ' + scl_file_wlz + \
              ' | WlzHistogramEqualiseObj -D | WlzGreyNormalise -u ' + \
              ' > ' + blr_file_wlz
        doCmd(cmd)
        cmd = 'WlzExtFFConvert -f wlz -F nii ' + \
                              ' -o ' + blr_file_nii + ' ' + blr_file_wlz
        doCmd(cmd)
      #}
    #}
  #}
  vrbMsg('End Section - prepare_data')
#}

if compute_registration: #{
  vrbMsg('Begin Section - compute registration')
  # Create composite objects in which the smoothed dapi and jag1 images
  # are added for all instances of each model.
  # The composite images are registered (using ANTs to determine a rigid 
  # body affine transform) to the first instance of the model.
  #
  # Create temp file paths
  n_tmp = 5
  tmp_wlz = []
  tmp_nii = []
  tmp_affine = tmp_file + '-crg-'
  for i in range(0,n_tmp): #{
    tmp_wlz.append(tmp_file + '-crg-' + str(i) + '.wlz')
    tmp_nii.append(tmp_file + '-crg-' + str(i) + '.nii')
  #}
  for mod in models: #{
    mod_idx = mod['index']
    vrbMsg('model index = ' + str(mod_idx))
    if mod_idx in model_idx_set: #{
      mod_asy_idx = mod['assays'][0]
      for asy_idx in mod['assays']: #{
        vrbMsg('assay index = ' + str(asy_idx))
        if asy_idx >= assay_idx_range[0] and \
           asy_idx <= assay_idx_range[1]: #{
          # Create images for registration. If the jag1 image channel
          # exists then use this to improve orientation
          asy = assays[asy_idx]
          nii_d_sb = blrFileChnNii(asy, Channel.dapi)
          wlz_d_sb = blrFileChnWlz(asy, Channel.dapi, True)
          wlz_j_sb = blrFileChnWlz(asy, Channel.jag1, True)
          if bool(wlz_j_sb): #{
            old_jag_dapi_combination = False
            if old_jag_dapi_combination: #{
              vrbMsg('affine registration using 4 x jagi + 1 x dapi')
              cmd = 'WlzGreySetRange -l 0 -u 255 -L 0 -U 49 ' + \
                    ' < ' + wlz_d_sb + ' > ' + tmp_wlz[1]
              doCmd(cmd)
              cmd = 'WlzGreySetRange -l 0 -u 255 -L 0 -U 199 ' + \
                    ' < ' + wlz_j_sb + ' > ' + tmp_wlz[2]
              doCmd(cmd)
              cmd = 'WlzImageArithmetic -a ' + tmp_wlz[1] + ' ' + tmp_wlz[2] + \
                    ' | ' + 'WlzGreyNormalise -u > ' + tmp_wlz[3]
              doCmd(cmd)
            else: #}{
              cmd = 'WlzBoundingBox ' + wlz_j_sb
              bb = doCmd(cmd).split()
              cmd = 'WlzThreshold -v 100 -H ' + wlz_j_sb + '|' \
                    'WlzGreyNormalise | ' + \
                    'WlzCutObjToBox -o ' + tmp_wlz[1] + \
                        ' -x ' + bb[0] + ',' + bb[3] + \
                        ' -y ' + bb[1] + ',' + bb[4] + \
                        ' -z ' + bb[2] + ',' + bb[5]
              doCmd(cmd)
              cmd = 'WlzImageArithmetic -a ' + tmp_wlz[1] + ' ' + wlz_d_sb + \
                    '| WlzGreyNormalise -u > ' + tmp_wlz[3]
              doCmd(cmd)
            #}
          else: #}{
            vrbMsg('affine registration using just dapi')
            cmd = 'WlzGreyNormalise -u ' + wlz_d_sb + ' > ' + tmp_wlz[3]
            doCmd(cmd)
          #}
          cmd = 'WlzExtFFConvert -f wlz -F nii -o ' + tmp_nii[3] + ' ' + \
                tmp_wlz[3]
          doCmd(cmd)
          if asy_idx == mod_asy_idx: #{
            # If this is the model assay then make it the registration target
            vrbMsg('model assay index = ' + str(mod_asy_idx))
            mod_nii_d_sb = nii_d_sb
            renameFile(tmp_wlz[3], tmp_wlz[0])
            renameFile(tmp_nii[3], tmp_nii[0])
          else: #}{
            # Check to see if a manual registration landmark file
            # exists, use it if it does otherwise compute an affine
            # registration using the blurred dapi/jag1 composite.
            lmk = landmarkFile(asy_idx, mod_asy_idx)
            if bool(lmk): #{
              vrbMsg('using landmark file ' + lmk)
              wlz_afn_tfm = wlzAfnTfmFile(mod, asy)
              cmd = 'WlzITKAffineTransform.py -O m -I l -o - ' + lmk + '|' + \
                    'WlzAffineTransformObj -3 -I -m - -N -T ' + ' ' + wlz_afn_tfm
              doCmd(cmd)
              cmd = 'WlzAffineTransformObj -3 -t ' + wlz_afn_tfm + \
                                         ' -o ' + tmp_wlz[4] + ' ' + wlz_d_sb
              doCmd(cmd)
              cmd = 'WlzExtFFConvert -o ' + tmp_nii[4] + ' ' + tmp_wlz[4]
              doCmd(cmd)
            else: #}{
              vrbMsg('computing rigid-affine registration.')
              # Check if manual registration landmarks exist
              cmd = 'ANTS 3 '
              if verbose: #{
                cmd = cmd + ' -v '
              #}
              cmd = cmd + \
                    ' -t Rigid ' + \
                    ' -o ' + tmp_affine + \
                    ' -m CC\\['+ tmp_nii[0] + ',' + tmp_nii[3] + ',1,4\\] ' + \
                    ' --rigid-affine true'
              doCmd(cmd)
            #}
            # Now register the scaled but not blurred dapi images using
            # an elastic transform with the composite rigid body transform
            # as the initial affine transform.
            vrbMsg('computing elastic registration.')
            cmd = 'ANTS 3 '   
            if verbose: #{    
              cmd = cmd + ' -v '
            #}
            out_tr = tfmFile(mod, asy)
            cmd = cmd + \
                  ' -o ' + out_tr + ' ' + \
                  ' -r Gauss\\[4,1\\] '
            if bool(lmk): #{
              cmd = cmd + \
                  ' -t Elast\\[1.0\\] ' + \
                    ' -m CC\\['+ mod_nii_d_sb + ',' + tmp_nii[4]+ ',1,4\\]'
            else: #}{
              cmd = cmd + \
                    ' -t Elast\\[0.4\\] ' + \
                    ' -m CC\\['+ mod_nii_d_sb + ',' + nii_d_sb + ',1,4\\] ' + \
                    ' -a ' + tmp_affine + 'Affine.txt'
            #}
            doCmd(cmd)
          #}
        #}
      #}
    #}
  #}
  vrbMsg('End Section - compute registration')
#}
  
# Apply the registration transforms to create registered assays
if register_models: #{
  vrbMsg('Begin Section - register model assays')
  n_tmp = 1
  tmp_wlz = []
  tmp_nii = []
  tmp_affine = tmp_file + '-crg-'
  for i in range(0,n_tmp): #{
    tmp_wlz.append(tmp_file + '-crg-' + str(i) + '.wlz')
    tmp_nii.append(tmp_file + '-crg-' + str(i) + '.nii')
  #}
  for mod in models: #{
    vrbMsg('model index = ' + str(mod['index']))
    for asy_idx in mod['assays']: #{
      vrbMsg('assay index = ' + str(asy_idx))
      mod_asy_idx = mod['assays'][0]
      if asy_idx >= assay_idx_range[0] and \
         asy_idx <= assay_idx_range[1]: #{
        asy = assays[asy_idx]
        tr_aff_file = re.sub('\.nii$', 'Affine.txt', tfmFile(mod, asy))
        tr_els_file = re.sub('\.nii$', 'Warp.nii',   tfmFile(mod, asy))
        if asy_idx == mod['assays'][0]: #{
          ref_asy = asy
          nii_ref_file = sclFileChnNii(asy, Channel.dapi)
        #}
        for img in asy['images']: #{
          chn = img['channel']
          vrbMsg('channel = ' + channel[chn])
          # register both the scaled and blurred scaled images
          for idx in range(0,2): #{
            vrbMsg('idx (0 scl, 1 blr) = ' + str(idx))
            if(idx == 0): #{
              nii_in_file = sclFileChnNii(asy, chn)
              wlz_in_file = sclFileChnWlz(asy, chn)
            else: #}{
              nii_in_file = blrFileChnNii(asy, chn)
              wlz_in_file = blrFileChnWlz(asy, chn, True)
            #}
            nii_out_file = regFileChnNii(asy, chn, idx > 0)
            wlz_out_file = regFileChnWlz(asy, chn, idx > 0)
            if asy_idx == mod_asy_idx: #{
              # If this is the first assay of the model just normalise the
              # image and wite it out.
              cmd = 'WlzGreyNormalise -u ' + wlz_in_file + ' > ' + wlz_out_file
              doCmd(cmd)
            else: #}{
              # If this is not the first assay of the model then transform
              # it using the model assay transform.
              wlz_afn_tfm = wlzAfnTfmFile(mod, asy)
              if not os.path.isfile(wlz_afn_tfm): #{
               wlz_afn_tfm = None
              #}
              if bool(wlz_afn_tfm): #{

                cmd = 'WlzAffineTransformObj -3 -t ' + wlz_afn_tfm + \
                                           ' -o ' + tmp_wlz[0] + ' ' + \
                                           wlz_in_file
                doCmd(cmd)
                cmd = 'WlzExtFFConvert -o ' + tmp_nii[0] + ' ' + tmp_wlz[0]
                doCmd(cmd)
                nii_in_file = tmp_nii[0]
              #}
              cmd = 'antsApplyTransforms '
              if verbose: #{    
                cmd = cmd + ' -v '
              #}
              cmd = cmd + \
                    ' -n Linear ' + \
                    ' -i ' + nii_in_file + ' ' + \
                    ' -r ' + nii_ref_file + ' ' + \
                    ' -o ' + nii_out_file + ' ' + \
                    ' -t ' + tr_els_file + ' ' + \
                    ' -t ' + tr_aff_file + ' '
              #}
              doCmd(cmd)
              cmd = 'WlzExtFFConvert -f nii -F wlz -o - ' + nii_out_file + \
                    ' | WlzGreyNormalise -u > ' + wlz_out_file
              doCmd(cmd)
            #}
          #}
        #}
      #}
    #}
  #}
  vrbMsg('End Section - register models')
#}

# Threshold assays, splitting the above threshold values into n_thr_band
# bands for each assay signal (ie non dapi) channel, then combine these to
# model summed band "occupancy" images
if create_gx_band_images: #{
  vrbMsg('Begin Section - create gene expression band images')
  tmp_files = []
  n_tmp_files = 3
  for i in range(0, n_tmp_files): #{
    tmp_files.append(tmp_file + '-thr-' + str(i) + '.wlz')
  #}
  for mod in models: #{ For each model
    mod_idx = mod['index']
    if mod_idx in model_idx_set: #{
      for asy_idx in mod['assays']: #{ For each assay of the model
        asy = assays[asy_idx]
        first = True
        for img in asy['images']: #{
          chn_idx = img['channel']
          if not chn_idx == Channel.dapi: #{
            thr_band_wid = 0.0
            thr_val_0 = 0
            thr_val = 0
            thr_val_0 = int(round(img['thresh_8bit']))
            thr_band_wid = (255.0 - thr_val_0) / n_thr_band
            vrbMsg('GX assay = ' + str(asy_idx) + \
                   ', channel = ' + channel[chn_idx] + \
                   ', thresh = ' + str(thr_val_0) + \
                   ', width = ' + str(thr_band_wid))
            # Create zero valued registered image
            cmd = 'WlzGreySetValue -g0 ' + \
                  regFileChnWlz(asy, Channel.dapi, False) + \
                  ' > ' + tmp_files[0] # t0
            doCmd(cmd)
            for thr_idx in range(0, n_thr_band): #{
              thr_val = int(round(thr_val_0 + thr_idx * thr_band_wid))
              cmd = 'WlzThreshold -H -v ' + str(thr_val) + ' ' + \
                     regFileChnWlz(asy, chn_idx, False) + ' | WlzDomain '
              if thr_band_dilation > 0: #{
                cmd = cmd + ' | WlzDilation -c26 -r' + str(thr_band_dilation)
              #}
              cmd = cmd + ' > ' + tmp_files[1] # WlzThreshold > t1
              doCmd(cmd)
              cmd = 'WlzGreyMask -m ' + str(thr_idx + 1) + ' ' + \
                    tmp_files[1] + ' ' + tmp_files[0] + \
                    ' > ' + tmp_files[2] # WlzGreyMask t1 t0 >t2
              doCmd(cmd)
              renameFile(tmp_files[2], tmp_files[0]) # mv t2 t0
            #}
            renameFile(tmp_files[0], thrFileAsyChnWlz(asy_idx, chn_idx)) # mv t0
          #}
        #}
      #}
    #}
  #}
#}

# Combine the assays into mean, max and normalised sum of histogram equalised
# assay expression models. Also combine the banded threshold images by
# adding them.
if combine_to_model: #{
  vrbMsg('Begin Section - combine to model')
  # Temporary files are indexed:
  #   0 zero filled cuboid image
  #   1 sum blurred
  #   2 max
  #   3 histogram equalised
  #   4 sum threshold bands
  #   5 sum threshold bands of blurred
  #   6 max threshold bands
  #   7 max threshold bands of blurred
  #   8 occupancy of lowest level thresholded images
  #   9 reserved
  #   others (11 12 13 14 15 16 17 18 19 20 21 22 23 24) are transient files 
  n_tmp_files = 25
  tmp_files = []
  for i in range(0, n_tmp_files): #{
    tmp_files.append(tmp_file + '-com-' + str(i) + '.wlz')
  #}
  for mod in models: #{
    # Create blank cuboid image into which to compose mean and max models
    mod_idx = mod['index']
    vrbMsg('Model index ' + str(mod_idx))
    if mod_idx in model_idx_set: #{
      vrbMsg('Creating zeroed image as basis for model ' + str(mod_idx))
      asy = assays[mod['assays'][0]]
      mod_sz_file = sclFileChnWlz(asy, Channel.dapi)
      cmd = 'WlzBoundingBox ' + mod_sz_file
      bb = doCmd(cmd).split()
      if dummy_run: #{
        bb = ['0', '1', '2', '3', '4', '5']
      #}
      bb_lim = ' -x ' + bb[0] + ',' + bb[3] + \
               ' -y ' + bb[1] + ',' + bb[4] + \
               ' -z ' + bb[2] + ',' + bb[5] + ' '
      cmd = 'WlzMakeRect -3  ' + bb_lim + \
            ' | WlzGreySetValue -g0 > ' + tmp_files[0]
      doCmd(cmd)
      # Now compose assay/channel images for model
      for chn in channel: #{
        chn_img_cnt = 0
        vrbMsg('Composing model images for channel ' + chn)
        chn_idx = matchChannel(chn)
        for idx in range(1, 10): #{
          copyFile(tmp_files[0], tmp_files[idx])
        #}
        for asy_idx in mod['assays']: #{
          asy = assays[asy_idx]
          for img in asy['images']: #{
            asy_chn_idx = img['channel']
            vrbMsg('Assay ' + str(asy_idx) + ' image channel ' + \
                   channel[asy_chn_idx])
            if asy_chn_idx == chn_idx: #{
              chn_img_cnt = chn_img_cnt + 1
              # Sum blurred images for channel
              vrbMsg('Computing sum of blurred images for channel ' + \
                     str(chn_idx))
              blr_img = regFileChnWlz(asy, chn_idx, True)
              vrbMsg('Using scaled blurred image ' + blr_img)
              cmd = 'WlzCutObjToBox ' + bb_lim + ' -o ' + tmp_files[11] + \
                    ' ' + blr_img
              doCmd(cmd)
              cmd = 'WlzImageArithmetic -a -o ' + tmp_files[12] + ' ' + \
                    tmp_files[1] + ' ' + tmp_files[11]
              doCmd(cmd)
              renameFile(tmp_files[12], tmp_files[1])
              # Max of scaled images for channel
              vrbMsg('Computing max of scaled images for channel ' + \
                     str(chn_idx))
              scl_img = regFileChnWlz(asy, chn_idx, False)
              vrbMsg('Using scaled image ' + scl_img)
              cmd = 'WlzCutObjToBox ' + bb_lim + ' -o ' + tmp_files[13] + \
                    ' ' + scl_img
              doCmd(cmd)
              cmd = 'WlzImageArithmetic -b 14 -o ' + tmp_files[14] + ' ' + \
                    tmp_files[2] + ' ' + tmp_files[13]
              doCmd(cmd)
              renameFile(tmp_files[14], tmp_files[2])
              # Sum of histogram equalised assays
              vrbMsg('Computing sum of histogram equalised blurred images ' + \
                     'for channel ' + str(chn_idx))
              blr_img = regFileChnWlz(asy, chn_idx, True)
              vrbMsg('Using scaled blurred image ' + blr_img)
              cmd = 'WlzThreshold -v20 -H ' + blr_img + ' | ' + \
                    'WlzHistogramEqualiseObj -D | ' + \
                    'WlzCutObjToBox ' + bb_lim + ' -o ' + tmp_files[15]
              doCmd(cmd)
              cmd = 'WlzImageArithmetic -a -o ' + tmp_files[16] + ' ' + \
                    tmp_files[3] + ' ' + tmp_files[15]
              doCmd(cmd)
              renameFile(tmp_files[16], tmp_files[3])
              # Sum of of gene expression thresholds
              if not chn_idx == Channel.dapi: #{
                vrbMsg('Computing sum of gene expression threshold images ' + \
                       'for channel ' + str(chn_idx))
                scl_img = thrFileAsyChnWlz(asy_idx, chn_idx)
                cmd = 'WlzImageArithmetic -a -o ' + tmp_files[17] + ' ' + \
                      scl_img + ' ' + tmp_files[4]
                doCmd(cmd)
                renameFile(tmp_files[17], tmp_files[4])
                vrbMsg('Computing sum of blurred gene expression ' + \
                       'threshold images for channel ' + str(chn_idx))
                scl_img = thrFileAsyChnWlz(asy_idx, chn_idx)
                cmd = 'WlzSepFilterObj -gd -m 2,2,2 -o ' + tmp_files[18] + \
                      ' ' + scl_img
                doCmd(cmd)
                cmd = 'WlzImageArithmetic -a -o ' + tmp_files[19] + ' ' + \
                      tmp_files[18] + ' ' + tmp_files[5]
                doCmd(cmd)
                renameFile(tmp_files[19], tmp_files[5])
                vrbMsg('Computing max of gene expression threshold images ' + \
                       'for channel ' + str(chn_idx))
                scl_img = thrFileAsyChnWlz(asy_idx, chn_idx)
                cmd = 'WlzImageArithmetic -b 14 -o ' + tmp_files[20] + ' ' + \
                      scl_img + ' ' + tmp_files[6]
                doCmd(cmd)
                renameFile(tmp_files[20], tmp_files[6])
                vrbMsg('Computing max of blurred gene expression ' + \
                       'threshold images for channel ' + str(chn_idx))
                scl_img = thrFileAsyChnWlz(asy_idx, chn_idx)
                cmd = 'WlzSepFilterObj -gd -m 2,2,2 -o ' + tmp_files[21] + \
                      ' ' + scl_img
                doCmd(cmd)
                cmd = 'WlzImageArithmetic -b 14 -o ' + tmp_files[22] + ' ' + \
                      tmp_files[21] + ' ' + tmp_files[7]
                doCmd(cmd)
                renameFile(tmp_files[22], tmp_files[7])
              #}
              # Occupancy of gene expression at/above threshold, this is the
              # occupancy of the first threshold band in the assay banded
              # threshold images
              if not chn_idx == Channel.dapi: #{
                vrbMsg('Computing occupancy of gene expression threshold ' + \
                       'images for channel ' + str(chn_idx))
                scl_img = thrFileAsyChnWlz(asy_idx, chn_idx)
                cmd = 'WlzThreshold -H -v1 ' + scl_img +  ' | ' + \
                      'WlzGreySetValue -g1 | ' + \
                      'WlzCutObjToBox ' + bb_lim + ' -o ' + tmp_files[23]
                doCmd(cmd)
                cmd = 'WlzImageArithmetic -a -o ' + tmp_files[24]  + ' ' + \
                      tmp_files[23] + ' ' + tmp_files[8]
                doCmd(cmd)
                renameFile(tmp_files[24], tmp_files[8])
              #}
            #}
          #}
        #}
        # Normalise the values and output
        if chn_img_cnt < 1: #{
          vrbMsg('No images in model ' + str(mod['index']) + \
                 ' in channel ' + channel[chn_idx])
        else: #}{
          # Normalise the blurred values and output
          vrbMsg('Normalising model images')
          mod_igm = modFileBlrChnWlz(mod['index'], chn_idx)
          cmd = 'WlzGreyNormalise -u ' + ' ' + tmp_files[1] + \
                ' | WlzThreshold -v20 -H  > ' + mod_igm
          doCmd(cmd)
          mod_igm = modFileSclChnWlz(mod['index'], chn_idx)
          cmd = 'WlzGreyNormalise -u ' + ' ' + tmp_files[2] + \
                ' | WlzThreshold -v20 -H > ' + mod_igm
          doCmd(cmd)
          mod_igm = modFileHeqChnWlz(mod['index'], chn_idx)
          cmd = 'WlzGreyNormalise -u ' + ' ' + tmp_files[3] + \
                ' | WlzThreshold -v20 -H > ' + mod_igm
          doCmd(cmd)
          if not chn_idx == Channel.dapi: #{
            mod_igm = modFileThrChnWlz(mod['index'], chn_idx, False, 0)
            renameFile(tmp_files[4], mod_igm)
            mod_igm = modFileThrChnWlz(mod['index'], chn_idx, True, 0)
            cmd = 'WlzConvertPix -t3 < ' + tmp_files[5] + ' > ' + mod_igm
            doCmd(cmd)
            mod_igm = modFileThrChnWlz(mod['index'], chn_idx, False, 1)
            renameFile(tmp_files[6], mod_igm)
            mod_igm = modFileThrChnWlz(mod['index'], chn_idx, True, 1)
            cmd = 'WlzConvertPix -t3 < ' + tmp_files[7] + ' > ' + mod_igm
            doCmd(cmd)
            mod_occ = modFileThrOccWlz(mod['index'], chn_idx, 1)
            renameFile(tmp_files[8], mod_occ)
          #}
        #}
      #}
    #}
  #}
  vrbMsg('End Section - combine to model')
#}

# Create point clouds for the models and registered images of the model
if create_point_clouds: #{
  vrbMsg('Begin Section - create point clouds')
  gam = 2.0
  pnt_cfg = '-D 1,1,1 -d 1 -G -g 64,255,' + str(gam) + ' '
  for mod in models: #{
    mod_idx = mod['index']
    vrbMsg('Model index ' + str(mod_idx))
    if mod_idx in model_idx_set: #{
      vrbMsg('Creating point clouds for model ' + str(mod_idx))
      for chn in channel: #{
        vrbMsg('Channel ' + chn)
        chn_idx = matchChannel(chn)
        if chn_idx > 0: #{
          in_file = modFileBlrChnWlz(mod['index'], chn_idx)
          if os.path.exists(in_file): #{
            pts_file = modPtsFileChnVtk(mod['index'], chn_idx)
            cmd = 'WlzPointsFromDomain ' + pnt_cfg + in_file +  \
                  ' | WlzExtFFConvert -f wlz -F vtk -o ' + pts_file + ' - '
            doCmd(cmd)
          else: #}{
            vrbMsg('Warning file ' + in_file + ' does not exist.')
          #}
        #}
      #}
      for asy_idx in mod['assays']: #{
        asy = assays[asy_idx]
        vrbMsg('Assay ' + str(asy_idx))
        for img in asy['images']: #{
          chn_idx = img['channel']
          vrbMsg('Channel ' + channel[chn_idx])
          if chn_idx == Channel.dapi: #{
            in_file = regFileChnWlz(asy, chn_idx, True)
            pts_file = regPtsFileChnVtk(asy, chn_idx)
            cmd = 'WlzPointsFromDomain ' + pnt_cfg + in_file + \
                  ' | WlzExtFFConvert -f wlz -F vtk -o ' + pts_file + ' - '
            doCmd(cmd)
          #}
        #}
      #}
    #}
  #}
  vrbMsg('End Section - create point clouds')
#}

# Create surfaces for the model dapi channel blurred images
if create_ref_surface: #{
  vrbMsg('Begin Section - create ref surfaces')
  tmp_dom_file = tmp_file + '-dom.wlz'
  for mod in models: #{
    mod_idx = mod['index']
    vrbMsg('Model index ' + str(mod_idx))
    if mod_idx in model_idx_set: #{
      vrbMsg('Model ' + str(mod['index']))
      in_file = modFileBlrChnWlz(mod['index'], Channel.dapi)
      srf_file = modSrfFileChnVtk(mod['index'], Channel.dapi)
      cmd = 'WlzThreshold -v 40 -H ' + in_file + \
            ' | WlzDomain | WlzDilation -c26 -r6 | WlzDomainFill ' + \
            ' | WlzErosion -c26 -r4 > ' + tmp_dom_file
      doCmd(cmd)
      cmd = 'WlzDomainToVTKSurf.py -o ' + srf_file + ' -f -m 5000 ' + \
            tmp_dom_file
      doCmd(cmd)
    #}
  #}
  vrbMsg('End Section - create ref surfaces')
#}

# Create a compound Woolz object for each of the models which has all the
# component (grey valued) channels
if create_compound: #{
  vrbMsg('Begin Section - create compound object')
  for mod in models: #{
    mod_idx = mod['index']
    vrbMsg('Model index ' + str(mod_idx))
    if mod_idx in model_idx_set: #{
      vrbMsg('Model ' + str(mod['index']))
      # Create blank object file for any missing channels
      blank_file = tmp_file + '-cpd-0.wlz'
      cmd = 'WlzGreySetValue -g0  ' + \
            modFileBlrChnWlz(mod['index'], Channel.dapi) + \
            ' > ' + blank_file
      doCmd(cmd)
      file_list = ''
      cpd_file = modCpdFile(mod['index'])
      for chn in channel: #{
        vrbMsg('Channel ' + chn)
        chn_idx = matchChannel(chn)
        if chn_idx > 0: #{
          in_file = modFileBlrChnWlz(mod['index'], chn_idx)
          if(not os.path.exists(in_file)): #{
            vrbMsg('Using blank in place of non-existing ' + in_file)
            in_file = blank_file
          #}
          file_list = file_list + ' ' + in_file
        #}
      #}
      cmd = 'cat ' + file_list + ' | WlzCompound -o ' + cpd_file
      doCmd(cmd)
    #}
  #}
  vrbMsg('End Section - create compound object')
#}

# Cut sections through all registered assay dapi channel blurred images
if cut_sections: #{
  vrbMsg('Begin Section - cut sections')
  for mod in models: #{
    mod_idx = mod['index']
    vrbMsg('Model index ' + str(mod_idx))
    if mod_idx in model_idx_set: #{
      first = True
      vrbMsg('Model ' + str(mod['index']))
      for asy_idx in mod['assays']: #{
        asy = assays[asy_idx]
        vrbMsg('Assay ' + str(asy['index']))
        in_file = regFileChnWlz(asy, Channel.dapi, True)
        if first: #{
          cmd = 'WlzCentreOfMass -b ' + in_file
          com = doCmd(cmd).split()
        #}
        for pln_idx in range(0, 3): #{
          pln_angles = {0: '0,0,0', 1: '90,0,0', 2: '90,90,0'}
          angles = pln_angles[pln_idx]
          sec_file = secFileChnTif(asy['index'], angles, Channel.dapi)
          cmd = 'Wlz3DGetSection -f ' + com[1] + ',' + com[2] + ',' + com[3] + \
                    ' -a ' + angles + ' ' + in_file + \
                    ' | WlzExtFFConvert -f wlz -F tif -o ' + sec_file + ' - '
          doCmd(cmd)
        #}
        first = False
      #}
    #}
  #}
  cmd = 'rm -f Working/sec-montage.png'
  doCmd(cmd)
  cmd = 'montage -label \'%t\' -tile 3x ' + \
                 '-size \'100%x100%\' Working/sec-* ' + \
                 'Working/sec-montage.png'
  doCmd(cmd)
  vrbMsg('End Section - cut sections')
#}

# Compute displacements
if compute_displacements: #{
  vrbMsg('Begin Section - compute displacements')
  for mod in models: #{
    mod_idx = mod['index']
    vrbMsg('Model index ' + str(mod_idx))
    if mod_idx in model_idx_set: #{
      all_dsp = ''
      vrbMsg('Model ' + str(mod['index']))
      for asy_idx in mod['assays']: #{
        if not asy_idx == mod['assays'][0]: #{
          asy = assays[asy_idx]
          vrbMsg('Assay ' + str(asy['index']))
          # Convert ANTs NIfTI transform to Woolz
          tfm = re.sub('\.nii$', 'Warp.nii',   tfmFile(mod, asy))
          dsp = dspFile(mod, asy)
          cmd = 'WlzExtFFConvert -o - -F wlz ' + tfm + ' | ' + \
                'WlzCompoundArrayToScalar -o ' + dsp
          doCmd(cmd)
          all_dsp = all_dsp + ' ' + dsp
        #}
      #}
      cmd = 'cat ' + all_dsp + ' | ' + \
            'WlzCompound -t 1 - | ' + \
            'WlzNObjsGreyStats -m -s - | ' + \
            'WlzExplode -p -b ' + modDspFileBody(mod_idx)
      doCmd(cmd)
    #}
  #}
#}

# Compute threshold band gene expression intersections
if compute_gx_intersect: #{
  vrbMsg('Begin Section - compute threshold band gene expression intersections')
  vrbMsg('Reading all threshold band gene expression images into array')
  errNum = w.enum__WlzErrorNum(w.WLZ_ERR_NONE)
  dapi_thr = 32
  occ_thr = 4                     # Use the lowest threshold band for occupancy
  hi = w.enum__WlzThresholdType(w.WLZ_THRESH_HIGH)
  # typ = 0: intersections of sum of threshold bands
  # typ = 1: intersections of max of threshold bands
  # typ = 2: intersections of occupancy values
  for typ in range(0,3): #{
    typ_str = ['sum of threshold bands', \
               'max of threshold bands',
               'occupancy values']
    vrbMsg('Computing intersections of ' + typ_str[typ])
    for mod in models: #{
      mod_idx = mod['index']
      vrbMsg('Model index ' + str(mod_idx))
      if mod_idx in model_idx_set: #{
        vrbMsg('Model ' + str(mod['index']))
        # Create blank object file for any missing channels
        blank_file = tmp_file + '-cgx-0.wlz'
        cmd = 'WlzGreySetValue -g0  ' + \
              modFileBlrChnWlz(mod_idx, Channel.dapi) + \
              ' > ' + blank_file
        doCmd(cmd)
        #
        thr_objs = [ None for i in range(len(channel))]
        isn_file_name = isnFile(mod_idx, typ)
        # HACK here for append rather than open
        # HACK isn_file = open(isn_file_name, 'a')
        isn_file = open(isn_file_name, 'w')
        print('channel_0, threshold_0, ' + \
              'channel_1, threshold_1, ' + \
              'vol_itersect, vol_union, jaccard', \
              file = isn_file)
        for chn_idx in range(0, len(channel)): #{
          if (chn_idx == Channel.unknown): #{
            pass
          else: #}{
            obj = None
            chn = channel[chn_idx]
            vrbMsg('Channel ' + chn)
            if (chn_idx == Channel.dapi): #{
              thr_file = modFileBlrChnWlz(mod_idx, chn_idx)
            else: #}{
              thr_file = modFileThrChnWlz(mod_idx, chn_idx, False, typ)
            #}
            vrbMsg('Reading ' + thr_file)
            if(not os.path.exists(thr_file)): #{
              vrbMsg('Warning ' + thr_file + ' does not exist, using blank file.')
              thr_file = blank_file
            #}
            try: #{
                fp = c.cast(libc.fopen(thr_file, 'rb'), c.POINTER(w.FILE))
                obj = w.WlzReadObj(fp, c.byref(errNum))
                libc.fclose(fp)
            except Exception as e: #}{
              vrbMsg('Failed to read Woolz object (errNum = ' + str(errNum) + ' exception = ' + str(e) + ')')
              pass
            #}
            if bool(obj):#{
              if (chn_idx == Channel.dapi): #{
                thr_objs[chn_idx] = w.WlzAssignObject(
                                    w.WlzThresholdI(obj, dapi_thr, hi,
                                                    c.byref(errNum)), None)
                w.WlzFreeObj(obj)
                obj = None
              else: #}{
                thr_objs[chn_idx] = w.WlzAssignObject(obj, None)
              #}
            else: #}{
              vrbMsg('File ' + thr_file + ' does not exist.')
              if(exit_on_cmd_error): #{
                sys.exit(1)
              #}
            #}
          #}
        #}
        vrbMsg('Computing intersections for all channel/threshold ' +
               'combinations')
        try: #{
          for chn_idx_0 in range(0, len(channel)): #{
            if (chn_idx_0 == Channel.unknown): #{
              pass
            else: #}{
              vrbMsg('Channel 0 ' + channel[chn_idx_0])
              max_0 = 1
              if not (chn_idx_0 == Channel.dapi): #{
                if (typ == 2): #{
                  max_0 = occ_thr
                else: #}{
                  max_0 = objMaxValue(thr_objs[chn_idx_0])
                #}
              #}
              for g_0 in range(1, max_0 + 1): #{
                obj_t_0 = None
                if (chn_idx_0 == Channel.dapi): #{
                  obj_t_0 = w.WlzAssignObject(thr_objs[chn_idx_0], None)
                else: #}{
                  obj_t_0 = w.WlzAssignObject(
                            w.WlzThresholdI(thr_objs[chn_idx_0],
                                            g_0, hi,
                                            c.byref(errNum)), None)
                #}
                name_0 = channel[chn_idx_0] + '_' + str(g_0)
                # HACK for partial range of channels (use with append)
                # for chn_idx_1 in range(0, Channel.mecom): #{
                # for chn_idx_1 in range(Channel.mecom, len(channel)): #{
                # HACK
                for chn_idx_1 in range(0, len(channel)): #{
                  if (chn_idx_1 == Channel.unknown): #{
                    pass
                  else: #}{
                    vrbMsg('Channel 1 ' + channel[chn_idx_1] + ' Channel 0 ' + channel[chn_idx_0])
                    max_1 = 1
                    if not (chn_idx_1 == Channel.dapi): #{
                      if (typ == 2): #{
                        max_1 = occ_thr
                      else: #}{
                        max_1 = objMaxValue(thr_objs[chn_idx_1])
                      #}
                    #}
                    for g_1 in range(1, max_1 + 1): #{
                      obj_t_1 = None
                      if (chn_idx_1 == Channel.dapi): #{
                        obj_t_1 = w.WlzAssignObject(thr_objs[chn_idx_1], None)
                      else: #}{
                        obj_t_1 = w.WlzAssignObject(
                                  w.WlzThresholdI(thr_objs[chn_idx_1],
                                                  g_1, hi,
                                                  c.byref(errNum)), None)
                      #}
                      if bool(errNum): #{
                        raise WlzError()
                      #}
                      name_1 = channel[chn_idx_1] + '_' + str(g_1)
                      # Compute volume of intersection.
                      obj_i = w.WlzAssignObject(
                              w.WlzIntersect2(obj_t_0, obj_t_1,
                                              c.byref(errNum)), None)
                      if bool(errNum): #{
                        raise WlzError()
                      #}
                      vol_i = w.WlzVolume(obj_i, c.byref(errNum))
                      if bool(errNum): #{
                        raise WlzError()
                      #}
                      # Compute volume of union.
                      obj_u = w.WlzAssignObject(
                              w.WlzUnion2(obj_t_0, obj_t_1,
                                          c.byref(errNum)), None)
                      if bool(errNum): #{
                        raise WlzError()
                      #}
                      vol_u = w.WlzVolume(obj_u, c.byref(errNum))
                      if bool(errNum): #{
                        raise WlzError()
                      #}
                      if vol_u < 1: #{
                        jaccard = 0.0
                      else: #}{
                        jaccard = float(vol_i) / float(vol_u)
                        if jaccard < 1.0e-06: #{
                          jaccard = 0.0
                        #}
                      #}
                      print(channel[chn_idx_0] + ', ' + 
                            str(g_0) + ', ' + 
                            channel[chn_idx_1] + ', ' +
                            str(g_1) + ', ' +
                            str(vol_i) + ', ' +
                            str(vol_u) + ', ' +
                            ('%1.6f'% jaccard),
                            file = isn_file)
                      w.WlzFreeObj(obj_i)
                      obj_i = None
                      w.WlzFreeObj(obj_u)
                      obj_u = None
                      w.WlzFreeObj(obj_t_1)
                      obj_t_1 = None
                    #}
                  #}
                #}
              #}
              w.WlzFreeObj(obj_t_0)
              obj_t_0 = None
            #}
          #}
        except: #}{
          err = w.WlzStringFromErrorNum(errNum, None)
          vrbMsg('Failed to compute intersections (' + err + ')')
          if exit_on_cmd_error: #{
            traceback.print_exc()
            sys.exit(1)
          #}
        #}
        isn_file.close()
        # Free Woolz objects
        for obj in thr_objs: #{
          if bool(obj): #{
            # TODO There may be a bug here as freeing this object
            # sometimes causes segfault
            # w.WlzFreeObj(obj)
            pass
          #}
        #}
      #}
    #}
  #}
#}

