import numpy as np
import pandas as pd 
import math
def getregressors(txtpath, workdir, TR, order, peak):

    # Open the txt file
    with open(txtpath, "r+") as f:
        d = f.readlines()

    # Remove the head infos 
    for i in range(11):
        d.pop(0)
    
    # Split values into colums
    dataset = []
    for str in d:
        row = str.split("\t")
        row.pop(3)
        dataset.append(row)
    
    # Give names to colums
    import pandas as pd 
    df = pd.DataFrame(dataset)
    df.columns = ['ppg','respiration','mri']
    t_scan, n_scan = pulse_detect(df['mri'])
    # Limit the length of t_scan by comparing the ends
    r_end = t_scan[n_scan-1] + 2*TR*10000
    while r_end <= len(df['mri']):
        try:
            df = df[0:r_end]
            break
        except ValueError:
            print("The physiological signal doesn't cover the whole scan")


    resp_phase = getphase_res(df, t_scan, peak, workdir)
    card_phase = getphase_ppg(df['ppg',workdir])

    # Calculate and Define the size of regressors
    m = order
    regressors = []
    for i in range(n_scan):
        t_reg = []
        for j in range(m):
            unit = [math.cos((j+1) * resp_phase[i]), math.sin((j+1) * resp_phase[i]), math.cos((j+1) * card_phase[i]), math.sin((j+1) * card_phase[i])]
            t_reg = np.concatenate((t_reg,unit))
        regressors.append(t_reg)
    regressors = pd.DataFrame(regressors)
    regressors.to_csv(workdir + '/regressors.csv', index = None)


def pulse_detect(mri):
    # Get the time stamp of each MRI volume by checking the pulse value
    # tag contains all the t_scan
    t_scan = mri
    pulse = t_scan[t_scan[:] > 3.1]
    index_p = pulse.index
    tag = []
    tag.append(index_p[0])
    ref = index_p[0]

    for i in range(len(index_p)):
        if i == len(index_p) - 1:
            break

        if index_p[i+1] -  ref > 50:
            tag.append(index_p[i+1])
            ref = index_p[i+1]
    n_scan = len(tag)   
    return t_scan, n_scan


def getphase_ppg(ppg, workdir):
    import heartpy as hp

    # Use heartpy to detect peaks
    working_data, measures = hp.process(ppg, 10000.0)

    #plot with different title
    fig = hp.plotter(working_data, measures, title='Heart Beat Detection on Noisy Signal')
    fig.savefig(workdir + '/heartbeat.png')

    peak = working_data.get("peaklist")

    return peak


def getphase_res(df, t_scan, peak, workdir):
    import numpy as np
    resp = np.zeros([len(df),2])
    resp[:,0] = df[:,1]

    # Mark the unreliable data segments
    resp[:,1] = np.where(resp[:,0] > 0.09, False, True)

    # Get the max and min value of respiration data
    max_value = np.max(resp[np.where(resp[:,1] == True),0])
    min_value = np.min(resp[np.where(resp[:,1] == True),0])
    range = max_value - min_value

    # Normalize the respiration data
    f = lambda x: (x-min_value)/(range)
    resp_norm = f(resp[np.where(resp[:,1] == True),0])
    hist = np.histogram(resp_norm, bins =100)

    num = len(hist[0])
    x = np.arange(num)
    hist[0]

    # Plot and Save the histogram of respiration
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(15, 5))
    plt.bar(x, height = hist[0], )
    fig.savefig(workdir + '/respiration_histogram.png')

    # Find the MRI pulse in the histogram and get the value
    hist = np.asarray(hist)

    # Calculate the value of each scan before regressors calculation
    val = []
    for i in t_scan:
        amp = resp[i,0]
        amp_norm = f(amp)
        val.append(phase_value(hist,amp_norm) * find_phase(peak,i[0]) * math.pi)
    return val


## Define the function that find the left amplitude
def find_left(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    if array[idx] > value:
        return idx-1
    else:
        return idx


## Define the function that calculate the value from histogram A/(A + B)
def phase_value(hist, amp_norm):
    idx = find_left(hist[1],amp_norm)
    num = len(hist[0])
    A = 0
    sum = 0
    for i in np.arange(num):
        if i <= idx:
            A = A + hist[0][i]
        sum = sum + hist[0][i]
    return A/sum


## Define the function that find the phase of respiration
## array here is resp_pulse and value here is mri pulse

# def find_leftphase(array, value):
#     array = np.asarray(array[:,0])
#     idx = (np.abs(array - value)).argmin()
#     if array[idx] > value:
#         return idx-1
#     else:
#         return idx
        

# Find the phase positive or negative with the help of 'peak'
# Peak is generated with bioPAC, containing the info of respiration state
def find_phase(array, value, peak):
    resp = np.asarray(array)
    sim_resp = resp[np.where(np.equal(resp[:,1],"Recovery") | np.equal(resp[:,1],"Expire Start"))]
    state = peak[find_left(resp[:,0],value/10000),1]

    if state == "Inspire Start" :
        state = peak[find_left(sim_resp[:,0],value/10000),1]
    if state == "Expire Start":
        return -1
    else:
        return 1