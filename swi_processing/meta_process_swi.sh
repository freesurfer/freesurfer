#!/bin/bash
# A script which calls process_swi for all the human swi scans

# MUSC Siemens
#./process_swi.sh siemens MUSC_Siemens ../0001_MUSC_20100825_MR/0001_MUSC_20100825_MR/SCANS/5-gre_swi_sag/DICOM/custom.INTR_MUSC_001.MR.INTR_DEV.5.1.20100825.141627.359000.lcgwu2.dcm ../0001_MUSC_20100825_MR/0001_MUSC_20100825_MR/SCANS/6-gre_swi_sag/DICOM/custom.INTR_MUSC_001.MR.INTR_DEV.6.1.20100825.141627.359000.19g1os2.dcm

# Cincinnati Philips
#./process_swi.sh philips CIN_Philips ../0001_CIN_08172010_MR/0001_CIN_08172010_MR/SCANS/501-swi/DICOM/custom.Cin_Philips_Dev_08_17_10_NormalViewDICOM_MRI_HUMAN_SUBJECT.MR.INTR_DEV.501.1.20100816.102456.1r5314z.dcm

# Dartmouth Philips
#./process_swi.sh philips DMS_Philips ../0001_DMS_20101102_MR/SCANS/601/DICOM/S601.00001

# UCSD GE
#./process_swi.sh ge UCSD_GE ../0001_UCSD_20101104_MR/SCANS/7/DICOM/i922826.MRDC.1

# BWH GE
./process_swi.sh ge BWH_GE ../0001_SWI/SCANS/11/DICOM/custom.0001_MR.MR.INTR_DEV.11.1.20101006.124443.ytl2lu.dcm

# DUKE GE
./process_swi.sh ge DUKE_GE ../0001_DUKE_20101027_MR1/SCANS/6/DICOM/custom.0001_MR.MR.INTR_DEV.6.1.20101027.131151.iruxgn.dcm

