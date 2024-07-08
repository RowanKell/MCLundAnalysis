import xlsxwriter
import uproot
import os

def convertData(xlsxFileName, rootFileName,multipleFiles):
    if(xlsxFileName == ""):
        xlsxFileName = 'xlsx/May23_500k.xlsx'
    if(rootFileName == ""):
        rootFileName = "for_box_500k.root"
    workbook = xlsxwriter.Workbook(xlsxFileName)
    worksheet = workbook.add_worksheet()
    kinematic_list = ["pT","Q2","x","z","R2_adjust","hadron","target"]
    kinematic_list_short = ["pT","Q2","x","z","R2_adjust"]
    # syntax: write(row,column,data)
    global_row_num = 0
    column = 1
    for kin in kinematic_list:
        worksheet.write(global_row_num,column,kin)
        column += 1;
    if(multipleFiles):
        file_dir = rootFileName #should be a directory
        num_files = len([name for name in os.listdir(file_dir) if not os.path.isdir(name)])
        file_names = [name for name in os.listdir(file_dir) if not os.path.isdir(name)]
    else:
        num_files = 1
        file_names = [rootFileName]
#     print(f"fileName: {file_names[0]}")

    tree_MC_list = []
    for i in range(num_files):
        try:
#             print(f"file_name: {file_names[i]}")
            tree_MC_list.append(uproot.open(file_dir + file_names[i]+ ":tree_MC"))
        except uproot.exceptions.KeyInFileError as e:
            print(f"exception: {e}\nexception for file {file_names[i]}; continuing")
            continue

    var_row = []
#     print(f"num_files: {num_files}")
    file_row_start = 1
    for i in range(num_files):
        j = 1
        for var in kinematic_list_short:
            print(f"file_row_start: {file_row_start}")
            row_num = file_row_start
            curr_array = tree_MC_list[i][var].array(library='np')
            for row in range(len(curr_array)):
                worksheet.write(row_num,j,curr_array[row])
                row_num += 1
            j += 1
        for row in range(len(tree_MC_list[0][kinematic_list_short[0]].array(library='np'))):
            worksheet.write(row_num,6,"pi+")
            worksheet.write(row_num,7,"proton")
        file_row_start += len(tree_MC_list[0][kinematic_list_short[0]].array(library='np'))
    workbook.close()
def binnedConvertData(xlsxFileName, rootFileName):
    if(xlsxFileName == ""):
        xlsxFileName = 'xlsx/binned_Lund_May22_with_R2.xlsx'
    if(rootFileName == ""):
        rootFileName = "file_1_all_bins_w_new_R2.root"
    workbook = xlsxwriter.Workbook(xlsxFileName)
    worksheet = workbook.add_worksheet()
    kinematic_list = ["x", "z", "Q2", "pT","R2","hadron","target"]
    # syntax: write(row,column,data)
    global_row_num = 0
    column = 1
    for kin in kinematic_list:
        worksheet.write(global_row_num,column,kin)
        column += 1;
#     file_dir = "/work/clas12/users/rojokell/MCLundAnalysis/OutputFiles/Files_Spring_24/May_22/"
    # num_files = len([name for name in os.listdir(file_dir) if not os.path.isdir(name)])
    # file_names = [name for name in os.listdir(file_dir) if not os.path.isdir(name)]
    num_files = 1
    file_names = [rootFileName]

    tree_x_list = []
    tree_z_h_list = [] 
    tree_qTdivQ_list = []
    for i in range(num_files):
        try:
            tree_x_list.append(uproot.open(file_names[i]+ ":tree_x_bins"))    
            tree_z_h_list.append(uproot.open(file_names[i]+ ":tree_z_h_bins"))    
            tree_qTdivQ_list.append(uproot.open(file_names[i]+ ":tree_qTQ_bins"))
        except uproot.exceptions.KeyInFileError as e:
            print(f"exception: {e}\nexception for file {file_names[i]}; continuing")
            continue

    xarray = np.array([np.array([np.zeros(7)] * 4)] * num_files)
    qTdivQarray = np.array([np.array([np.zeros(9)] * 5)] * num_files)
    zarray = np.array([np.array([np.zeros(7)] * 4)] * num_files)

    xkinematics = np.array(["z_h", "Q2", "pT","R2"])
    zkinematics = np.array(["x", "Q2", "pT","R2"])
    qTdivQkinematics = np.array(["x", "z_h", "Q2", "pT","R2"])

    xbins = np.array([0.1,0.13,0.16,0.19,0.235,0.3,0.5])
    zbins = np.array([0.35,0.43,0.49,0.55,0.62,0.7,0.83])
    qTdivQbins = np.array([0.1,0.3,0.5,0.8,1.1,1.5,2,2.5,3])
    global_row_num += 1 #increment past labels
    #column numbers based on this list: ["pT","Q2","x","z","R2","hadron","target"]
    z_col_num = 4
    pT_col_num = 1
    Q2_col_num = 2
    x_col_num = 3

    for i in range(num_files):
        z_iter = 0
        for var in zkinematics:
            zarray[i][z_iter] = tree_z_h_list[i][var].array(library='np')
            z_iter += 1
        x_iter = 0
        for var in xkinematics:
            xarray[i][x_iter] = tree_x_list[i][var].array(library='np')
            x_iter += 1
        qTdivQ_iter = 0
        for var in qTdivQkinematics:
            qTdivQarray[i][qTdivQ_iter] = tree_qTdivQ_list[i][var].array(library='np')
            qTdivQ_iter += 1

    xarray_t = np.array([np.array([np.zeros(4)] * 7)] * num_files)
    qTdivQarray_t = np.array([np.array([np.zeros(5)] * 9)] * num_files)
    zarray_t = np.array([np.array([np.zeros(4)] * 7)] * num_files)
    for i in range(num_files):
        xarray_t[i] = np.transpose(xarray[i])
        zarray_t[i] = np.transpose(zarray[i])
        qTdivQarray_t[i] = np.transpose(qTdivQarray[i])
    z_pos = [1,3,4,5]
    row_num = 1
    for i in range(num_files):
        for bin_num in range(len(xarray_t[i])):
            for var_num in range(len(xarray_t[i][bin_num])):
                worksheet.write(row_num,1,xbins[bin_num]) #write the x value from bin
                worksheet.write(row_num,var_num+2,xarray_t[i][bin_num][var_num]) #write the other kin. from LundAnalysis binning for this bin
            row_num += 1 #inc the row after each bin
        for bin_num in range(len(zarray_t[i])):
            for var_num in range(len(zarray_t[i][bin_num])):
                worksheet.write(row_num,2,zbins[bin_num]) #write the x value from bin
                worksheet.write(row_num,z_pos[var_num],zarray_t[i][bin_num][var_num]) #write the other kin. from LundAnalysis binning for this bin
            row_num += 1 #inc the row after each bin
        for bin_num in range(len(qTdivQarray_t[i])):
            for var_num in range(len(qTdivQarray_t[i][bin_num])):
                try:
                    worksheet.write(row_num,var_num+1,qTdivQarray_t[i][bin_num][var_num]) #write the other kin. from LundAnalysis binning for this bin
                except:
                    if math.isnan(qTdivQarray_t[i][bin_num][var_num]):
                        worksheet.write(row_num,var_num+1,0)
            row_num += 1 #inc the row after each bin
        for row in range(row_num - 1):
            worksheet.write(row+1,6,"pi+")
            worksheet.write(row+1,7,"proton")
    workbook.close()
    
