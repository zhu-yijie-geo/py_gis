# -*- coding: utf-8 -*-
"""
Created on Tue May 16 20:14:05 2023

@author: Yijie Zhu
"""
import random,string,os,time,datetime,shutil
import numpy as np
import pandas as pd
from datetime import datetime
import pandas as pd
import time
import arcpy
arcpy.env.overwriteOutput = True
from arcpy import env
from arcpy.sa import *



basedir="D:/AAAG_516"
###You may change the basedir to your prefered location
"""
This script gives example to perform 3 runs. Each run will take about 20 mins depending on your PC configuration.
The code here does not employ multiprocessing. The excution time will rely on the speed of the single core.  
Before we start, make sure in your master folder: 
    1. put ibtracs data in the folder 'ibtracs'
    2. put boundary/shp data in the folder 'data/Boundary' 
"""
#####-pre_process
###-EXTRACT 1900-2020/6-h
print('PROCESS STARTS')
st00=time.time()

geo_db='raster'
##create geodatabase to save raster files
arcpy.management.CreateFileGDB(basedir, geo_db, "CURRENT")
arcpy.env.workspace ='{}/{}'.format(basedir,'ibtracs')
ibtracs='{}/{}/{}'.format(basedir,'ibtracs','IBTrACS_original.shp')
##create a folder to save output tracks 
if not os.path.exists('{}/{}'.format(basedir,'track')):
    os.makedirs('{}/{}'.format(basedir,'track'))
outdir='{}/{}'.format(basedir,'track')
outname='P6.shp'
print('Filtering TC tracks')
##from ibtracs, select only 6 hour track points and limit time to 1900-2020
arcpy.analysis.Select(ibtracs, '{}/{}'.format(outdir,outname), "hour IN (0, 6, 12, 18) And min = 0 And SEASON >= 1900 And SEASON <= 2020")

### Project to NAD1983-17N
workspace="{}/track/".format(basedir)
arcpy.management.Project(workspace+"P6.shp", workspace+"P6_17.shp", 
                          'PROJCS["NAD_1983_UTM_Zone_17N",GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-81.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]', "NAD_1983_To_WGS_1984_1", 'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]', "NO_PRESERVE_SHAPE", None, "NO_VERTICAL")



"""
Now starts random run. The example here shows a set of random run for 20km, 50km, and a reference run 
Create random TC locations within a distance from the besttrack.
It will take around 20 mins to run each case depending on your PC configuration. 
"""
basefolder=basedir

run_list=['r20','r50','rc0']
run_list_r=['20000','50000','0.1']
for i, r in enumerate(run_list):
    st0=time.time()
    print('Start the {} set run'.format(r))
    arcpy.analysis.PairwiseBuffer(workspace+'P6_17', workspace+'P6_17{}'.format(r), "{} Meters".format(run_list_r[i]), "NONE", None, "PLANAR", "0 Meters")
    
    arcpy.management.CreateRandomPoints(workspace, "P6_17{}_rdm".format(r), workspace+"P6_17{}".format(r), "0 0 250 250", 1, "0 Meters", "POINT", 0)
    ##random points generated does not have storm info. Therefore we need to join storm info to the new file
    arcpy.management.JoinField(workspace+"P6_17{}_rdm.shp".format(r), "CID", workspace+"P6_17{}.shp".format(r), "FID", None)
    sc=arcpy.da.SearchCursor(workspace+"P6_17{}_rdm.shp".format(r), "SID")
    SIDLIST=[]
    for x in sc:
        if x not in SIDLIST:
            SIDLIST.append(x)
        else: continue
    infile=workspace+"P6_17{}_rdm.shp".format(r)
    out=workspace
    ##link points to lines
    print('-Linking TC Points to Lines')
    outname="L6_17{}_rdm.shp".format(r)
    arcpy.CreateFeatureclass_management(out,outname,"POLYLINE",None,"DISABLED","DISABLED",'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]];-400 -400 1111948722.22222;-100000 10000;-100000 10000;8.98315284119521E-09;0.001;0.001;IsHighPrecision', '', 0, 0, 0, '')
    
    arcpy.AddField_management(out+outname,"p6_ID_i","DOUBLE","","","10")
    arcpy.AddField_management(out+outname,"p6_ID_t","DOUBLE","","","10")
    arcpy.AddField_management(out+outname,"sid","TEXT","","","50")
    arcpy.AddField_management(out+outname,"Length","DOUBLE","","","3")
    arcpy.AddField_management(out+outname,"Vi","DOUBLE","","","3")
    arcpy.AddField_management(out+outname,"Vt","DOUBLE","","","3")
    arcpy.AddField_management(out+outname,"Pi","DOUBLE","","","3")
    arcpy.AddField_management(out+outname,"Pt","DOUBLE","","","3")
    arcpy.AddField_management(out+outname,"TIME_i","TEXT","","","50")
    arcpy.AddField_management(out+outname,"TIME_t","TEXT","","","50")  
    myflds=["Shape@XY","SID","USA_WIND","USA_PRES","ISO_TIME",'FID']
    sc1=arcpy.da.SearchCursor(infile, myflds)
    ic=arcpy.da.InsertCursor(out+outname,["Shape@",'p6_ID_i','p6_ID_t',"sid","Length","Vi","Vt","Pi","Pt","TIME_i","TIME_t"])
    
    ALLLIST=[]
    i=0
    for x in sc1:
        ALLLIST.append(x)
    del sc1
    for x in ALLLIST:
        if i<len(ALLLIST)-1:
            if x[1]==ALLLIST[i+1][1]:
                P1=arcpy.Point(x[0][0],x[0][1])
                P2=arcpy.Point(ALLLIST[i+1][0][0],ALLLIST[i+1][0][1])
                array=arcpy.Array([P1,P2])
                polyline = arcpy.Polyline(array)
                Length=polyline.length
                sid=x[1]
                Vi=x[2]
                Vt=ALLLIST[i+1][2]
                Pi=x[3]
                Pt=ALLLIST[i+1][3]
                TIME_i=x[-2]
                TIME_t=ALLLIST[i+1][-2]
                p6_ID_i=x[-1]
                p6_ID_t=ALLLIST[i+1][-1]
                ic.insertRow([polyline,p6_ID_i,p6_ID_t,sid,Length,Vi,Vt,Pi,Pt,TIME_i,TIME_t])
                i=i+1
            else:
                i=i+1
                continue
    del ic
    
    arcpy.management.DefineProjection(workspace+"L6_17{}_rdm.shp".format(r), 'PROJCS["NAD_1983_UTM_Zone_17N",GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-81.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]')
    ######Generate_Points_Along_Lines
    st=time.time()
    infile=workspace+"L6_17{}_rdm.shp".format(r)
    outfile=workspace+"L6_17{}_rdm_p1.shp".format(r)
    arcpy.management.GeneratePointsAlongLines(infile, outfile, "PERCENTAGE", None, 1, None)
    print('-Finished generating 100 segments along {}'.format(infile))
    print('-Output file saved to {}'.format(outfile))
    
    ###Intersect with Florida
    area="{}/data/Boundary/".format(basefolder)
    track='{}/track/'.format(basefolder)
    arcpy.management.Project(area+"FL.shp", area+"FL17.shp", 
                             'PROJCS["NAD_1983_UTM_Zone_17N",GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-81.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]', "NAD_1983_To_WGS_1984_1", 'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]', "NO_PRESERVE_SHAPE", None, "NO_VERTICAL")
    
    arcpy.management.Project(area+"US.shp", area+"US17.shp", 
                             'PROJCS["NAD_1983_UTM_Zone_17N",GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-81.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]', "NAD_1983_To_WGS_1984_1", 'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]', "NO_PRESERVE_SHAPE", None, "NO_VERTICAL")
    
     
    infiles=[track+"L6_17{}_rdm_p1.shp".format(r),area+"FL17.shp"]
    outfile=track+"L6_17{}_rdm_p1_FL.shp".format(r)
    arcpy.analysis.Intersect(infiles, outfile, "ALL", None, "INPUT")
    #### append storm info to L6_17{}_rdm_p1_FL.shp
    arcpy.management.JoinField(outfile, "ORIG_FID", track+"L6_17{}_rdm.shp".format(r), "FID", None)
    
    ####6_h point Intersect with US###
    ##This is for filtering different scenarios later
    
    p6='{}/track/P6_17{}_rdm.shp'.format(basefolder,r)
    infiles=[p6,area+"US.shp"]
    outfile=track+"P6_17{}_rdm_US.shp".format(r)
    arcpy.analysis.Intersect(infiles, outfile, "ALL", None, "INPUT")
    
    et=time.time()
    elapsed_time = (et - st0)/60 
    print('---Execution time for preprocessing run {}:'.format(r), elapsed_time, 'minutes')
    
    
    ####Export to xls for future use or other research purposes
    ####
    print('-Exporting preprocessed results to spreadsheets')
    st=time.time()
    export_list=['L6_17{}_rdm.shp'.format(r),"L6_17{}_rdm_p1_FL.shp".format(r), "L6_17{}_rdm_p1.shp".format(r), "P6_17{}_rdm_US.shp".format(r)]
    folder=basefolder
    for j in export_list:
        arcpy.conversion.TableToExcel(track+j, '{}/{}.xlsx'.format(folder,j[:-4]), "NAME", "CODE")
    et=time.time()
    elapsed_time = (et - st)/60 
    print('---Execution time for exporting preprocessed results to spreadsheets: ', elapsed_time, 'minutes')
    

    

    
    
   
      
    """
    this module calculates inland wind based on scenarios
    """
    print('-Start calculating inland wind')
    st=time.time()
    folder=basefolder
    
    df_L6p1FL = pd.read_excel('{}/L6_17{}_rdm_p1_FL.xlsx'.format(basefolder,r))
    
    df_P6US = pd.read_excel('{}/P6_17{}_rdm_US.xlsx'.format(basefolder,r))
    
    df=df_L6p1FL
    ###sort
    df.sort_values("FID_L6_17r",inplace=True)
     
    ##calculete the sequence (based on where the point locates in the 100 segments)
    def sequence(row):  
        return (row['FID_L6_17r']%99)+1 
    df["sequence"] = df.apply(lambda row: sequence(row), axis=1)
    ###
    ###create time for each track point based on the sequence
    df['time_i']=pd.to_datetime(df['TIME_i'])
    
    df["time_p1"] = pd.to_datetime(df['TIME_i']) + pd.to_timedelta(df['sequence']*216, unit='s')
    
    #####CHECK IF position i/t inland
    df_P6US.sort_values("FID_P6_17r",inplace=True)
    ####  i
    df = pd.merge(
        left=df,
        right=df_P6US['FID_P6_17r'],
        left_on='p6_ID_i',
        right_on='FID_P6_17r',
        how='left'
    )
    df['i_inland'] = df['FID_P6_17r']
    df.drop('FID_P6_17r', inplace=True, axis=1)
    df['i_inland'] = df['i_inland'].fillna(0)
    
    ####  t
    
    df = pd.merge(
        left=df,
        right=df_P6US['FID_P6_17r'],
        left_on='p6_ID_t',
        right_on='FID_P6_17r',
        how='left'
    )
    df['t_inland'] = df['FID_P6_17r']
    df.drop('FID_P6_17r', inplace=True, axis=1)
    df['t_inland'] = df['t_inland'].fillna(0)
    
    df.loc[df['t_inland'] >0]
    
    ####### FOUR SCENARIOS
    def scenario(row):
        if row['i_inland'] >0 and row['t_inland'] >0:
            return 1 ###in_in
        elif row['i_inland'] ==0 and row['t_inland'] >0:
            return 2 ###out_in
        elif row['i_inland'] >0 and row['t_inland'] ==0:
            return 3 ###in_out
        else:
            return 4 ###out_out
    
    df['scenario'] = df.apply(scenario, axis=1)
    
    ########
    #drop if v=0
    vt0_idx=df.loc[df['Vt']==0].index
    df.drop(vt0_idx, inplace=True)
    ##
    
    
    ##calculate decay constant
    df['a']=np.log(df['Vi']/df['Vt'])/6
    ##

    ####### inland v based on decay constant and on FOUR SCENARIOS
    def inland_v(row):
        if row['scenario'] ==1 or row['scenario'] ==2:
            return row['Vi']* math.e ** (-row['a']*(((row['time_p1']-row['time_i']).total_seconds())/3600))#exponential
        elif row['scenario'] ==3 and row['a'] >=0:
            return row['Vi'] 
        elif row['scenario'] ==3 and row['a'] <0:
            return row['Vi']* math.e ** (-row['a']*(((row['time_p1']-row['time_i']).total_seconds())/3600))#exponential
        elif row['scenario'] ==4:
            return row['Vi']-((row['Vi']-row['Vt'])/6)*((((row['time_p1']-row['time_i']).total_seconds())/3600))#linear interpolation
        else:
            return 999
    df['inland_v']=df.apply(inland_v, axis=1)

    
    #drop where v<34
    vt34_idx=df.loc[df['inland_v']<34].index
    df.drop(vt34_idx, inplace=True)
    
    
    ### Calculate ACE relative to 64 kt
    df['inland_v2']=df['inland_v']**2
    df['ACE']=df['inland_v2']/((64**2)*100)
    
    ##export to a csv file to be joined back to GIS shapefiles
    df.to_csv('{}/cal_out{}.csv'.format(folder,r))
    et=time.time()
    elapsed_time = (et - st)/60
    
    print('---Execution time for calculateing inland wind and exporting to spreadsheets: ', elapsed_time, 'minutes')
    run_time = (et - st00)/60
    print('Run Time Since Start: ',run_time,' minutes')
    ###create pivot table to summarize TC ACE per county per year
    print('-Start calculating County return period')
    st=time.time()
    df['Time_Year']=df['time_p1'].dt.year
    table = pd.pivot_table(df,index=['Time_Year'],columns=['County'],values=['ACE'],aggfunc=np.sum, fill_value=0)
    
    county_list=[]
    for i, col in enumerate(table.columns):
        county_list.append(col[1])
    
    
    ###any occurences
    ###here we have 121 years of data
    count_ace64=table[table['ACE']>=0.01].count()['ACE']
    ace64_return=list(121/count_ace64)
    
    count_ace96=table[table['ACE']>=0.0225].count()['ACE']
    ace96_return=list(121/count_ace96)
    
    
    # 1h of occurences
    count_ace641h=table[table['ACE']>=0.167].count()['ACE']
    ace64_return1h=list(121/count_ace641h)
    
    count_ace961h=table[table['ACE']>=0.375].count()['ACE']
    ace96_return1h=list(121/count_ace961h)
    ######
    
    ###create pivot table to summarize TC v_max per county per year
    table = pd.pivot_table(df,index=['Time_Year'],columns=['County'],values=['inland_v'],aggfunc=np.max, fill_value=0)
    count_v64=table[table['inland_v']>=64].count()['inland_v']
    v64_return=list(121/count_v64)
    count_v96=table[table['inland_v']>=96].count()['inland_v']
    v96_return=list(121/count_v96)
    
    
    ##merge into a dataframe and export
    ##output file can be joined back to the FL shapefile for plotting or conducting further analysis
    df_county=pd.DataFrame({'County': county_list, 'ace64_return': ace64_return, 'ace96_return': ace96_return,
                            'ace64_return1h':ace64_return1h, 'ace96_return1h':ace96_return1h,
                            'v64_return': v64_return, 'v96_return':v96_return})
    df_county.replace([np.inf, -np.inf], 999, inplace=True)
    
    df_county.to_csv('{}/county_out{}.csv'.format(basefolder, r))
    et=time.time()
    elapsed_time = (et - st)/60 
    print('-Finished Analysis Part. Execution time {}:'.format(r), elapsed_time, 'minutes')
    run_time = (et - st00)/60
    print('Run Time Since Start: ',run_time,' minutes')

    
    print('-Start generating raster files for plotting moving sum and moving average')


    infile="/track/L6_17{}_rdm_p1_FL.shp".format(r)
    outfile_34="/track/L6_17{}_rdm_p1_FLcal.shp".format(r)
    area="{}/data/Boundary/".format(basedir)
    arcpy.analysis.Select(folder+infile, folder+outfile_34, "Vt <> 0")
    arcpy.management.JoinField(folder+outfile_34, "FID_L6_17r", "{}/cal_out{}.csv".format(folder,r), "FID_L6_17r", 
                                "a;scenario;inland_v;ACE")
    

    
    ####moving average for v
    outfolder="{}/{}.gdb/".format(basedir,geo_db)
    mv34="mv_34_1050{}".format(r)
    ###output resolution 10km, search radius 50km x 50km
    out_raster = arcpy.sa.PointStatistics(folder+outfile_34, "inland_v", 10000, "Rectangle 50000 50000 MAP", "MEAN"); 
    out_raster.save(outfolder+mv34)
    ####resampled to 1km
    mv34_r='mv_34_1050_r1{}'.format(r)
    arcpy.management.Resample(outfolder+mv34, outfolder+mv34_r, "1000 1000", "BILINEAR")
    ####extract to only include Fl
    mv34_r_FL='mv34_1050_r1_FL{}'.format(r)
    out_raster = arcpy.sa.ExtractByMask(outfolder+mv34_r, area+"FL.shp"); 
    out_raster.save(outfolder+mv34_r_FL)
    
    ####moving SUM for ACE

    ms34="ms_34_1050{}".format(r)
    ###output resolution 10km, search radius 50km x 50km
    out_raster = arcpy.sa.PointStatistics(folder+outfile_34, "ACE", 10000, "Rectangle 50000 50000 MAP", "SUM"); 
    out_raster.save(outfolder+ms34)
    ####resampled to 1km
    ms34_r='ms_34_1050_r1{}'.format(r)
    arcpy.management.Resample(outfolder+ms34, outfolder+ms34_r, "1000 1000", "BILINEAR")
    ####extract to only include Fl
    ms34_r_FL='ms34_1050_r1_FL{}'.format(r)
    out_raster = arcpy.sa.ExtractByMask(outfolder+ms34_r, area+"FL.shp"); 
    out_raster.save(outfolder+ms34_r_FL)
    et=time.time()
    elapsed_time = (et - st0)/60
    print('-Finished generating raster files for plotting moving sum and moving average')
    print('Finished the {} set run; Total Execution time: '.format(r), elapsed_time, ' minutes')
    run_time = (et - st00)/60
    print('Run Time Since Start: ',run_time,' minutes')
    
et00=time.time()   
elapsed_time = (et00 - st00)/60 
print('CONGRATULATIONS! Whole Process Finished. Total Processing Time: ', elapsed_time, ' minutes')


