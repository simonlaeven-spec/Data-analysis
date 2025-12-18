import numpy as np
import matplotlib.pyplot as plt
from  scipy.integrate import simpson
from scipy.integrate import trapezoid
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from tkinter import filedialog
import os
from scipy.optimize import curve_fit
import h5py



#to make exe: pyinstaller --onefile --windowed .\delay_scan.py     in map stage (cd python workspace\stage)


def on_closing():
    root.destroy()  # Destroys the window
    root.quit()     # Stops the mainloop


root = tk.Tk()
root.title("Delay Scan")
root.geometry("800x600+100+100")



folderlength = 0



root.protocol("WM_DELETE_WINDOW", on_closing) #necessary to stop the program when closing the window


def instructions():
    window = tk.Toplevel(root)
    window.title("Instructions")
    window.geometry("400x300+150+150")
    instructions_text = ("- Select H5 file with scanning data.\n"
                         "- Inspect mass spectra by entering file number and clicking 'plot file'.\n"
                         "- Choose target peaks and fitting tolerance. Fits between target ± tolerance fit. \n"
                            "- Click 'Fit gaussian' to fit the peaks. A plot of the currently selected file pops up.\n"
                            "- Check whether fitted peaks are correct and check integration bounds. \n"
                            "- Plot times integrates between the integration bounds numerically for all files and plots areas vs time delay.\n"
                            "- Branching ratio plots the fraction of each target peak over the total integrated counts vs time delay.\n"


    )
    label = tk.Label(window, text=instructions_text, wraplength=380, justify="left")
    label.pack(pady=10, padx=10)




def browse_h5_file():
    '''Lets user browse h5 file'''
    
    file_path = filedialog.askopenfilename(filetypes=[("H5 files", "*.h5")])
    
    return file_path


def h5_file():
    '''reads out h5 file, returns arrays with for each scan the masses, magnitudes, keys and parameter values'''
    global file_loaded_label
    global root
    masses = []
    magnitudes = []
    filepath = browse_h5_file()
    
    file_loaded_label.config(text = "Loading file...")
    root.update_idletasks()
    with h5py.File(filepath, 'r') as f:
        #print(list(f.keys()))

        rawdata = f['Rawdat']
        parameters = f['Parameters'][()]
        print(parameters)

        scan_start, scan_stop, scan_step = np.abs(parameters.item()[:3])
         
        
        print((scan_stop-scan_start)/scan_step)
        parameter_values = np.linspace(scan_start, scan_stop, int((scan_stop-scan_start)/scan_step)+1)
        

        traces = {}

        for key in rawdata.keys():
            if key.startswith('P0'):
                group = rawdata[key]

                if 'Trace' in group:    
                    traces[key] = group['Trace'][:]
        
        #print(traces)
        i = 0
        for key, trace in traces.items():
            i += 1
            trace = np.array(trace).reshape(-1)
            #print(trace.shape)
            #plot_trace(trace,key)

            sampling_frequency = 32000000 #Hz
            
            magnitude_pos = np.absolute(np.fft.rfft(trace, len(trace))[1:]) #amplitudes corresponding to frequencies          
            frequencies_pos = np.fft.rfftfreq(len(trace), d=1/sampling_frequency)[1:] #frequency spectrum
            
            
            
            
            
            #plot_freq(frequencies_pos, magnitude_pos,key)

            B = 7 #T
            
            ms_hfre = 1.079e+8  
            ms_beta = -1.857e+2
            
            '''first line is theoretical way of determining mass array. second line introduces ms_beta for a better fit. B*1.535611*10**7 = 1.075e8'''
            #mz = (B*1.535611*10**7)/frequencies_pos
            mz = (ms_hfre/(frequencies_pos-ms_beta))

            mz_0_500 = mz[mz <= 500] #array with values going from high to low
            print (mz_0_500)
            magnitude_pos_0_500 = magnitude_pos[-len(mz_0_500):]

            mz_0_500 = mz_0_500[::-1] #array ordered from low to high (necessary for integration)
            magnitude_pos_0_500 = magnitude_pos_0_500[::-1] #array ordered from low to high
            masses.append(mz_0_500)
            magnitudes.append(magnitude_pos_0_500)

            #plot_ms(mz_0_500, magnitude_pos_0_500,key)

    return masses, magnitudes, rawdata.keys(), parameter_values





def plot_file():
    '''plots the mass spectrum of the selected file in the tkinter window'''
    global file_nr
    global all_masses
    global all_counts
    global parameter_values
    i = int(file_nr.get()) -1
    parameter_values_ms = np.array(parameter_values)*1e6
    fig, ax = plt.subplots()
    ax.plot(all_masses[i], all_counts[i], label = f"Time delay: {parameter_values_ms[i]} " r"($\mu$s)")
    ax.set_title(f"Mass spectrum for file {i+1}")
    ax.set_xlabel( 'mass (amu)')
    ax.set_ylabel('counts')
    ax.legend(loc = 'upper right')
    #create a figure and axis for the plot
    fig.tight_layout()
    canvas = FigureCanvasTkAgg(fig, master=root)  
    canvas.draw()
    canvas_widget = canvas.get_tk_widget()
    canvas_widget.place(x=325, y = 180, width=480, height=350)
    #embed the plot in the tkinter window

    toolbar_frame = tk.Frame(root)
    toolbar_frame.place(x = 325, y = 550, width= 480,height=30 )  # just above the plottoolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
    toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
    toolbar.update() 
    #add a toolbar for the plot
    





    
def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))


def fit_gaussian(x, y, target_peak, fit_tolerance):
    '''fits a gaussian to the peaks around target_peak with bounds defined by fit_tolerance. Retures a list of [amplitude, center, width], except when fit fails, then returns None'''
    fit_bounds = (target_peak - fit_tolerance, target_peak + fit_tolerance)
    # Select data in the fit range
    mask = (x >= fit_bounds[0]) & (x <= fit_bounds[1]) #fit bounds is a tuple (lower, upper)
    x_fit = x[mask]
    y_fit = y[mask]

    # Initial guess: amplitude, center, width
    a0 = np.max(y_fit) #max value in y_fit
    x0 = x_fit[np.argmax(y_fit)] #x value where y is max
    sigma0 = (fit_bounds[1] - fit_bounds[0]) / 5

    try: 
        popt, pcov = curve_fit(gaussian, x_fit, y_fit, p0=[a0, x0, sigma0], maxfev=50000)
    except RuntimeError:
        return None
    a_fit, x0_fit, sigma_fit = popt

    return [a_fit, x0_fit, sigma_fit]




def fit_func():
    '''fits gaussians to the three target peaks. Takes total counts, by summing over all files. '''
    global target1
    global target2
    global target3
    global tolerance
    global folderlength
    
    all_counts_total = np.zeros(len(all_counts[0])) 
    for i in range(folderlength): #sums all counts of each file.
        all_counts_total = all_counts_total + all_counts[i]
    fit_params = [] #list of fit parameters for each peak
    fit_bounds = [] #bounds of region where trying to fit
    fit_succesful = [] #keeps track of which peaks were succesfully fitted
    fit_tolerance = float(tolerance.get()) #tolerance for fitting, for determining fit bounds

    '''Following part performs fits to determine the actual peak positions'''
    for i in range(3):
        target_peak = float([target1.get(), target2.get(), target3.get()][i])#peak to fit
        fit_bounds.append(target_peak - float(fit_tolerance))#fit bounds array, 6 elements. left and right bound for each of the 3 peaks
        fit_bounds.append(target_peak + float(fit_tolerance))

        fit_result = fit_gaussian(all_masses[0], all_counts_total, float(target_peak), 0.1) #fit gaussian to the target peak
        if fit_result is None:
            a_fit, x0_fit, sigma_fit = 0.0, target_peak, 0 #keep center at target peak if fit fails
            fit_succesful.append(False)#list that keeps track of which peaks were succesfully fitted. list of 3x3 elements
        else:
            a_fit, x0_fit, sigma_fit = fit_result
            fit_succesful.append(True)

        if x0_fit >= target_peak + fit_tolerance or x0_fit <= target_peak - fit_tolerance:
            #if the fitted center is outside the fit bounds, consider the fit failed
            a_fit, x0_fit, sigma_fit = 0.0, target_peak, 0
            fit_succesful[-1] = False  #update last element to False

        fit_params.append([a_fit, x0_fit, sigma_fit])#list of 3 lists, each with fit parameters for each peak 
    fit_params = np.array(fit_params)
    return fit_params, fit_succesful, fit_bounds








def background(fit_params, fit_succesful):
    '''determines background levels for each of the three target peaks based on fitted peak positions.
    Determines background by taking mean of counts in regions defined by tolerance_bg and outside integration bounds defined by tolerance_int'''
    global tolerance_bg
    global tolerance_int
    global folderlength
    global all_masses
    global all_counts
    
    backgrounds_list = []
    mask_bg = []
    x_coord_bg = [] #list containing x coordinates of background regions for plotting. list of 3 arrays with coordinates
    bg_bounds = [] #list containing bounds for background of all 3 peaks. format: [left peak one, right peak one, left peak 2 etc...]
    int_bounds = [] #list containing integration bounds for all 3 peaks

    for i in range(3):
        
        '''determining the bg bounds for peak i'''
        x0_fit = fit_params[i][1]  #center of the fitted peak. i is index of peak.

        lower_bound = x0_fit - float(tolerance_int.get()) #fit_params[i][1] is the center of the fit peak for target i. sets bounnds for integration. 
        upper_bound = x0_fit + float(tolerance_int.get()) #are used in determining the background outside of the integration part

        bg_lower = x0_fit - float(tolerance_bg.get())        #background bounds
        bg_upper = x0_fit + float(tolerance_bg.get())

        mask_bg = ((all_masses[0] >= bg_lower) & (all_masses[0] < lower_bound)) | ((all_masses[0] > upper_bound) & (all_masses[0] <= bg_upper)) #mask for background region
        x_coord_bg.append(all_masses[0][mask_bg])
        bg_bounds.append(bg_lower)
        bg_bounds.append(bg_upper)
        int_bounds.append(lower_bound)
        int_bounds.append(upper_bound)
        ''''''
        backgrounds_peak = []
        for index in range(folderlength):
            if not fit_succesful[i]:
                backgrounds_peak.append(0.0)  #if fit failed, set background to 0
                continue
        


            y_background_region = all_counts[index][mask_bg] #array of counts in the background region
            background_level = np.mean(y_background_region) #mean of counts in bg region
        
            backgrounds_peak.append(background_level)# list with 3 arrays. first array contains all bg values for first peak.
        backgrounds_list.append(backgrounds_peak)
    return backgrounds_list, x_coord_bg, bg_bounds, int_bounds






def plot_fit():
    '''plots the fit results for the selected file'''
    global file_nr
    global all_masses
    global all_counts
    global fit_params
    global fit_bounds
    global backgrounds_list

    

    masses = all_masses[int(file_nr.get()) - 1]
    counts = all_counts[int(file_nr.get()) - 1]


    #for plotting:
    
 #for background plotting
    backgrounds_list_for_plot_0 = np.linspace(backgrounds_list[0][int(file_nr.get()) - 1], backgrounds_list[0][int(file_nr.get()) - 1], len(x_coord_bg[0])) #horizontal line at background level
    backgrounds_list_for_plot_1 = np.linspace(backgrounds_list[1][int(file_nr.get()) - 1], backgrounds_list[1][int(file_nr.get()) - 1], len(x_coord_bg[1]))
    backgrounds_list_for_plot_2 = np.linspace(backgrounds_list[2][int(file_nr.get()) - 1], backgrounds_list[2][int(file_nr.get()) - 1], len(x_coord_bg[2]))  
    

    plt.close('all')  #close previous plots
    plt.figure(figsize=(8, 5))
    
    plt.plot(masses, counts, label="Data")
    plt.plot(fit_bounds, [0,0,0,0,0,0], 'ro', label="Fit bounds")

    plt.plot(x_coord_bg[0], backgrounds_list_for_plot_0, 'k--')
    plt.plot(x_coord_bg[1], backgrounds_list_for_plot_1, 'k--')
    plt.plot(x_coord_bg[2], backgrounds_list_for_plot_2, 'k--')

    plt.plot(bg_bounds, [0,0,0,0,0,0], 'go', label="Background bounds")
    plt.plot(int_bounds, [0,0,0,0,0,0], 'mo', label="Integration bounds")

    if fit_succesful[0]:
        plt.axvline(x=fit_params[0][1], color='orange', linestyle='--', label=f'Peak1 = {fit_params[0][1]} ')
    if fit_succesful[1]:
        plt.axvline(x=fit_params[1][1], color='purple', linestyle='--', label=f'Peak2 = {fit_params[1][1]} ')
    if fit_succesful[2]:    
        plt.axvline(x=fit_params[2][1], color='brown', linestyle='--', label=f'Peak3 = {fit_params[2][1]} ') 



    plt.xlabel("m/z")
    plt.ylabel("Intensity")
    plt.title(f"Gaussian Fits to Target Peaks file {int(file_nr.get())}")
    plt.legend(loc='upper right')
    plt.show(block = False)











def integrate(index, peak_number):
    global fit_params
    global backgrounds_list
    global fit_succesful
    global target1
    global target2
    global target3

    masses = all_masses[index]
    counts = all_counts[index]
    
    lower_bound = fit_params[peak_number -1][1] - float(tolerance_int.get()) 
    upper_bound = fit_params[peak_number -1][1] + float(tolerance_int.get())

    bg = backgrounds_list[peak_number - 1][index]#background level for the target peak. first part returns array of 3 elements, second part selects correct background value
    counts_corrected_target = counts - bg #background subtraction
    counts_corrected_target[counts_corrected_target < 0] = 0 #sets negative values to 0 after background subtraction



    integrating_between = (masses >= lower_bound) & (masses <= upper_bound) #returns array with true/false with true when masses is between lower and upper bound
    
    x_peak = masses[integrating_between]
    y_peak = counts_corrected_target[integrating_between]
    
    area1 = simpson(y_peak, x_peak) 
    area2 = trapezoid(y_peak, x_peak) #2 different methods to calculate area under curve
    area = (area1+area2)/2
    return area

def integrate_button_func(index):
    global target1
    global target2
    global target3
    global fit_succesful
    
    if fit_succesful[0]:    
        area1 = integrate(index, 1)
    else:
        area1 = 0.0
    if fit_succesful[1]:
        area2 = integrate(index, 2)
    else:
        area2 = 0.0
    if fit_succesful[2]:
        area3 = integrate(index, 3)
    else:
        area3 = 0.0


    return area1, area2, area3



def timeplot():
    
    global folderlength
    global all_masses
    global all_counts
    global name1
    global name2
    global name3
    global parameter_values

    areas = []
    for i in range(folderlength): 
        masses = all_masses[i]
        counts = all_counts[i]
        area1, area2, area3 = integrate_button_func(i)
        areas.append([area1, area2, area3])
    areas = np.array(areas)
    parameter_values_ms = np.array(parameter_values)*1e6
    plt.close('all')  #close previous plots
    plt.figure(figsize=(8, 5))
    plt.plot(parameter_values_ms, areas[:,0], label=f"{name1.get()}")
    plt.plot(parameter_values_ms, areas[:,1], label=f"{name2.get()}")
    plt.plot(parameter_values_ms, areas[:,2], label=f"{name3.get()}")
    plt.scatter(parameter_values_ms, areas[:,0], s = 4)
    plt.scatter(parameter_values_ms, areas[:,1], s = 4)   
    plt.scatter(parameter_values_ms, areas[:,2], s = 4)
    plt.grid()
    plt.xlabel(u"Time delay (\u03bcs)")
    plt.ylabel("Integrated area (a.u.)")
    plt.title("Integrated peak areas vs time delay")
    plt.legend(loc='upper right')
    plt.show(block = False)

def branching_ratio():
    global folderlength
    global all_masses
    global all_counts
    global parameter_values
    global name1
    global name2
    global name3

    ratios1 = []
    ratios2 = []
    ratios3 = []
    for i in range(folderlength): 
        masses = all_masses[i]
        counts = all_counts[i]
        area1, area2, area3 = integrate_button_func(i)
        try:
            ratio3 = area3/(area1 + area2 + area3)
        except ZeroDivisionError:
            ratio3 = 0.0

        ratios3.append(ratio3)

        try:
            ratio1 = area1/(area1 + area2 + area3)
        except ZeroDivisionError:
            ratio1 = 0.0
        ratios1.append(ratio1)

        try:
            ratio2 = area2/(area1 + area2 + area3)
        except ZeroDivisionError:
            ratio2 = 0.0
        ratios2.append(ratio2)

    parameter_values_ms = np.array(parameter_values)*1e6
    plt.figure(figsize=(8, 5))
    plt.plot(parameter_values_ms, ratios1, label=f"Fraction {name1.get()} / ({name1.get()} + {name2.get()} + {name3.get()})")
    plt.scatter(parameter_values_ms, ratios1, s = 4)
    plt.plot(parameter_values_ms, ratios2, label=f"Fraction {name2.get()} / ({name1.get()} + {name2.get()} + {name3.get()})")
    plt.scatter(parameter_values_ms, ratios2, s = 4)
    plt.plot(parameter_values_ms, ratios3, label=f"Fraction {name3.get()}  / ({name1.get()} + {name2.get()} + {name3.get()})")
    plt.scatter(parameter_values_ms, ratios3, s = 4)
    plt.grid()
    plt.xlabel(u"Time delay (\u03bcs)")
    plt.ylabel("Branching Ratio")
    plt.title("Branching Ratio vs time delay")
    plt.legend(loc='upper right')
    plt.show(block = False)



def fit_button_clicked():
    global fit_params
    global fit_succesful
    global fit_bounds
    global backgrounds_list
    global x_coord_bg
    global bg_bounds
    global int_bounds

    fit_params, fit_succesful, fit_bounds = fit_func()
    backgrounds_list, x_coord_bg, bg_bounds, int_bounds = background(fit_params, fit_succesful)
    plot_fit()
    





file_nr_label = ttk.Label(root, text="File number:")
file_nr_label.place(x=15, y=100, width=100, height= 20)

file_nr = tk.StringVar(value="1")  #default value
textbox_fn = tk.Entry(root, textvariable=file_nr)
textbox_fn.place(x = 130, y = 100, width=100, height=20)

plot_file_button = tk.Button(root, text="Plot file", command=plot_file)
plot_file_button.place(x=250, y=100, width=100, height=20)


instructions_button = ttk.Button(root, text="Manual", command= instructions)
instructions_button.place(x=15, y=20, width=150  , height=30)

target1_label = ttk.Label(root, text="Target 1:")
target1_label.place(x=15, y=130, width=100, height= 20)

target1 = tk.StringVar(value="181")  #default value
textbox_t1 = ttk.Entry(root, textvariable=target1)
textbox_t1.place(x = 130, y = 130, width=100, height=20)

name1 = tk.StringVar(value="Ta")  #default value
textbox_name1 = ttk.Entry(root, textvariable=name1)
textbox_name1.place(x = 240, y = 130, width=50, height=20)


file_loaded_label = ttk.Label(root,text = "File not loaded" )
file_loaded_label.place(x = 15, y = 70, width = 100, height=20)



target2_label = ttk.Label(root, text="Target 2:")
target2_label.place(x=15, y=160, width=100, height= 20)

target2 = tk.StringVar(value = "197")
textbox_t2 = ttk.Entry(root, textvariable=target2)
textbox_t2.place(x = 130, y = 160, width=100, height=20)

name2 = tk.StringVar(value="TaO")  #default value
textbox_name2 = ttk.Entry(root, textvariable=name2)
textbox_name2.place(x = 240, y = 160, width=50, height=20)


target3_label = ttk.Label(root, text="Target 3:")
target3_label.place(x=15, y=190, width=100, height= 20)

target3 = tk.StringVar(value = "212.9")
textbox_t3 = ttk.Entry(root, textvariable=target3)
textbox_t3.place(x = 130, y = 190, width=100, height=20)

name3 = tk.StringVar(value=r"TaO$_2$")  #default value
textbox_name3 = ttk.Entry(root, textvariable=name3)
textbox_name3.place(x = 240, y = 190, width=50, height=20)

tolerance_label = ttk.Label(root, text="Tolerance fit:")
tolerance_label.place(x=15, y=220, width=100, height= 20)

tolerance = tk.StringVar(value = "0.3")
textbox_tol = ttk.Entry(root, textvariable=tolerance)
textbox_tol.place(x = 130, y = 220, width=100, height=20)


tolerance_integrating_label = ttk.Label(root, text="Tol intregation:")
tolerance_integrating_label.place(x=15, y=250, width=100, height= 20)

tolerance_int = tk.StringVar(value = "0.05")
textbox_tol_int = ttk.Entry(root, textvariable=tolerance_int)
textbox_tol_int.place(x = 130, y = 250, width=100, height=20)



tolerance_bg_label = ttk.Label(root, text="Tol background:")
tolerance_bg_label.place(x=15, y=280, width=100, height= 20)

tolerance_bg = tk.StringVar(value = "0.5")
textbox_tol_bg = ttk.Entry(root, textvariable=tolerance_bg)
textbox_tol_bg.place(x = 130, y = 280, width=100, height=20)






fit_button = ttk.Button(root, text="Fit gaussian", command= fit_button_clicked)
fit_button.place(x=15, y=310, width=150  , height=30)


#integrate labels + button:
        
branching_ratio_button = ttk.Button(root, text="Plot branching ratio", command= branching_ratio)
branching_ratio_button.place(x=15, y=340, width=150  , height= 30)

time_button = ttk.Button(root, text="Plot times", command= timeplot)
time_button.place(x=15, y=370, width=150  , height= 30)





all_masses, all_counts, all_keys, parameter_values = h5_file()
folderlength = len(parameter_values)
fit_button_clicked()
plot_file()
file_loaded_label.config(text = "File loaded!")
root.update_idletasks()



root.mainloop()