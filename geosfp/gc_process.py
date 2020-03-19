"""
Created on March 18, 2020

@author: Yi Wang
"""

import datetime
import os


#
#------------------------------------------------------------------------------
#
def replace_line(ori_list, search_key, new_line):
    """
    """

    if new_line[-1] != '\n':
        new_line = new_line + '\n'

    for i in range(len(ori_list)):
        line = ori_list[i]
        if search_key in line:
            ori_list[i] = new_line

#
#------------------------------------------------------------------------------
#
def create_job(yyyymmdd, root_run_dir, geosfp_dir, 
        root_output_dict):
    """ Create a job for prcessing GEOS-FP data
    (ywang, 03/18/2020)

    Parameters
    ----------
    yyyymmdd : str
        Date, '20180501' for example.
    root_run_dir : str
        Parent directory of run directory.
    geosfp_dir : str
        Raw GEOS-FP directory.
    root_output_dict : dict
        Keys: Output type, '==> 2 x 2.5 output' for example.
        Values: Parent directory of output

    """

    if root_run_dir[-1] != '/':
        root_run_dir = root_run_dir + '/'

    if geosfp_dir[-1] != '/':
        geosfp_dir = geosfp_dir + '/'

    # run directory
    run_dir = root_run_dir + yyyymmdd
    if os.path.isdir(run_dir):
        os.system('rm -rf ' + run_dir)
    os.system('mkdir -p ' + run_dir)

    # directory of current script
    file_dir = os.path.dirname(os.path.abspath(__file__))

    # read gc_process_job.sh
    f_run_in = open(file_dir + '/gc_process_job.sh', 'r')
    job_lines = f_run_in.readlines()
    f_run_in.close()

    # 'time ./GeosFpDriver0.x < date_filename' in gc_process_job.sh
    replace_line(job_lines, 'time ./GeosFpDriver0.x < date_filename', 
            'time ./GeosFpDriver0.x < ' + yyyymmdd)

    # 'cd work_directory' in gc_process_job.sh
    replace_line(job_lines, 'cd work_directory',
            'cd ' + run_dir)

    # output gc_process_job.sh
    f_run_out = open(run_dir + '/gc_process_job.sh', 'w')
    for line in job_lines:
        f_run_out.write(line)
    f_run_out.close()
    os.system('chmod +x ' + run_dir + '/gc_process_job.sh')

    # yyyymmdd file
    f = open(run_dir + '/' + yyyymmdd, 'w')
    f.write(yyyymmdd)
    f.close()

    # read GeosFpDriver.input
    fin = open(file_dir + '/GeosFpDriver.input', 'r')
    gfdi = fin.readlines()
    fin.close()

    # geosfp_dir
    for i in range(len(gfdi)):
        if ('==> Local Raw Data Path' in gfdi[i]):
            gfdi[i+1] = geosfp_dir + '\n'

    # root_output_dict
    for out_type in root_output_dict:

        for i in range(len(gfdi)):

            if (out_type in gfdi[i]):

                # output dir
                out_dir = root_output_dict[out_type]
                if out_dir[-1] != '/':
                    out_dir = out_dir + '/'
                out_dir = out_dir + 'YYYY/MM/'

                # logical variable
                gfdi[i+1] = 'T\n'

                # output dir
                gfdi[i+3] = out_dir + '\n'
                gfdi[i+4] = out_dir + '\n'

                # create out_dir_real
                out_dir_real = out_dir.split('/')
                out_dir_real[-2] = yyyymmdd[4:6]
                out_dir_real[-3] = yyyymmdd[0:4]
                out_dir_real = '/'.join(out_dir_real)
                if not os.path.isdir(out_dir_real):
                    os.system('mkdir -p ' + out_dir_real)


    # output GeosFpDriver.input
    fout = open(run_dir + '/GeosFpDriver.input', 'w')
    for line in gfdi:
        fout.write(line)
    fout.close()

    # copy files
    file_list = ['GeosFpDriver0.x', 'GeosFpTemplateFile.nc', 
            'weights_025x03125_to_05x0625.txt', 
            'weights_025x03125_to_2x25.txt',
            'weights_025x03125_to_4x5.txt']
    for filename in file_list:
        os.system('cp ' + file_dir + '/' + filename 
                + ' ' + run_dir)

#
#------------------------------------------------------------------------------
#
def process_geosfp_month(yyyymm, root_run_dir, root_geosfp_dir, 
        root_output_dict):
    """ Process GEOS-FP one-month asm tavg1 data.
    (ywang, 03/19/2020)
    yyyymm : str
        Month, '201805' for example.
    root_run_dir : str
        Parent directory of run directory.
    root_geosfp_dir : str
        Parent directory of raw GEOS-FP directory.
    root_output_dict : dict
        Keys: Output type, '==> 2 x 2.5 output' for example.
        Values: Parent directory of output

    """

    if root_run_dir[-1] != '/':
        root_run_dir = root_run_dir + '/'

    if root_geosfp_dir[-1] != '/':
        root_geosfp_dir = root_geosfp_dir + '/'

    currDate_D = datetime.datetime.strptime(yyyymm + '01', '%Y%m%d')

    #cwd_dir = os.getcwd()

    while True:

        # Date
        currDate = str(currDate_D)
        yyyy = currDate[0:4]
        mm   = currDate[5:7]
        dd   = currDate[8:10]
        yyyymmdd = yyyy + mm + dd
        print('process ' + yyyymmdd)

        # run directory
        run_dir = root_run_dir + yyyymmdd

        # geosfp directory
        geosfp_dir = root_geosfp_dir + 'Y' + yyyy + '/M' + mm + '/D' + dd

        # create job
        create_job(yyyymmdd, root_run_dir, geosfp_dir, root_output_dict)

        # submit job 
        os.chdir(run_dir)
        os.system('qsub gc_process_job.sh')

        # Next day
        nextDate_D = currDate_D + datetime.timedelta(days=1)
        if (str(nextDate_D)[5:7] != mm):
            break
        currDate_D = nextDate_D












#
#------------------------------------------------------------------------------
#
