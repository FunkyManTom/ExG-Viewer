import os
import urllib.request
import io
import csv
import uproot
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import matplotlib.patches as ptc
from matplotlib.widgets import Slider,RadioButtons,Button
import numpy as np
from skimage.transform import resize

root = tk.Tk()
root.withdraw()
currdir = os.getcwd()

root_mat=uproot.open(filedialog.askopenfilename(filetypes=[("Root file",".root")],parent=root, initialdir=currdir, title='Select your ExG Matrix'))
print("\nIndex\tKey")
for idx,key in enumerate(root_mat.keys()):
    print("\n"+str(idx)+":\t"+key)
root_name=root_mat.keys()[int(input("\nChoose key index: "))]
primary_vals_rs=root_mat[root_name].to_numpy()[0].transpose()#np.load(filedialog.askopenfilename(filetypes=[("Numpy Matrix",".npy")],parent=root, initialdir=currdir, title='Select your (square!) ExG Matrix'))

if len(primary_vals_rs)!=len(primary_vals_rs[0]):
    print("!! NON SQUARE MATRIX DETECTED, REDUCING TO AXIS OF LOWER RESOLUTION !!")
    lower_dim=min(len(primary_vals_rs),len(primary_vals_rs[0]))
    primary_vals_rs = resize(primary_vals_rs, (lower_dim, lower_dim), anti_aliasing=False)

primary_Ex=root_mat[root_name].to_numpy()[1]
if len(primary_Ex)>len(root_mat[root_name].to_numpy()[2]):
    primary_Ex=root_mat[root_name].to_numpy()[2]

def f(os1,os2):
    #print(len(primary_vals_rs),len(primary_vals_rs[0]))
    diag_arr = np.diagonal(primary_vals_rs, offset=-np.min([os1, os2]))
    arr_max = len(diag_arr)

    for i in range(np.min([os1, os2]), np.max([os1, os2])):
        tmp_arr = np.diagonal(primary_vals_rs, offset=-i)
        diag_arr = diag_arr + np.pad(tmp_arr, (0, arr_max - len(tmp_arr)))

    return diag_arr

size_exg=len(primary_Ex)

reader_levels="init"
jp_pref="init"
iso_yn=input("\nLoad isotope levels (y/n)?: ").lower()
if iso_yn=="y" or iso_yn=="yes" or iso_yn=="ye":
    isotope=input("\nEnter isotope: ")
    isotope=isotope.lower()
    digits = ''.join(ch for ch in isotope if ch.isdigit())
    letters = ''.join(ch for ch in isotope if ch.isalpha())
    isotope=digits+letters

    url="https://www-nds.iaea.org/relnsd/v1/data?fields=levels&nuclides="+isotope
    req=urllib.request.Request(url)
    req.add_header('User-Agent','Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:77.0) Gecko/20100101 Firefox/77.0')

    with urllib.request.urlopen(req) as response:
        data_levels = response.read().decode('utf-8')

    csv_file = io.StringIO(data_levels)
    reader_levels = csv.DictReader(csv_file)
    jp_pref=input("\nEnter spin parity (e.g. 2+, leave blank for all levels): ")

plt.ion()
plt.style.use('dark_background')
# Create initial data
initial_var1 = 1
initial_var2 = 0
x = primary_Ex[0:len(f(initial_var1,initial_var2))]
y = f(initial_var1,initial_var2)

x_mat=[-initial_var1,size_exg-initial_var1]
y_mat=[0,size_exg]

# Set up the figure and the main plot axis
fig, (ax1,ax2) = plt.subplots(figsize=(15,8),ncols=2)

primary_imshow=ax1.imshow(np.log(primary_vals_rs+1),origin="lower",cmap="viridis")
if iso_yn=="y" or iso_yn=="yes" or iso_yn=="ye":
    for row in reader_levels:
        level_e=row["energy"]
        level_jp=row["jp"]
        if jp_pref=="":
            tmp_val = np.where(primary_Ex == min(primary_Ex, key=lambda k: abs(k - float(level_e))))[0][0]
            ax1.plot([-tmp_val, size_exg - tmp_val], [0, size_exg],"--w",label=level_e)
        elif level_jp==jp_pref:
            tmp_val=np.where(primary_Ex==min(primary_Ex,key=lambda k:abs(k-float(level_e))))[0][0]
            ax1.plot([-tmp_val,size_exg-tmp_val],[0,size_exg],"--w",label=level_e)

diag1, =ax1.plot(x_mat,y_mat,c="red")
diag2, =ax1.plot(x_mat,y_mat,c="red")
d_area= ax1.add_patch(ptc.Polygon([(0,-initial_var1),(0,-initial_var2),(size_exg,size_exg-initial_var1),(size_exg,size_exg-initial_var2)],color="r",alpha=0.3))
ax1.set_ylim(0,size_exg-5)
ax1.set_xlim(0,size_exg-5)
tdat1=ax1.text(size_exg-size_exg/4,0.04*size_exg,s=str(int(np.round(primary_Ex[initial_var1])))+" keV",size=24,c="white")
tdat2=ax1.text(size_exg-size_exg/4,0.1*size_exg,s=str(int(np.round(primary_Ex[initial_var2])))+" keV",size=24,c="white")

#ax2.set_yscale("log")
ax2.set_xlim(0,np.round(primary_Ex[size_exg-initial_var1]))
ax2.set_ylim(5,np.max(f(initial_var1,initial_var2))*1.1)
ax2.yaxis.set_label_text("Counts",rotation=0)
ax2.yaxis.set_label_coords(-.05,1)
ax2.xaxis.set_label_text("Eg [keV]")
ax2.xaxis.set_label_coords(.95,-.05)
# Adjust the bottom to make room for the slider
plt.subplots_adjust(bottom=0.25)
line, = ax2.plot(x, y,c="r")

# Create the slider's axis and the Slider object itself
slider_ax1 = plt.axes([0.13, 0.15, 0.34, 0.03])  # [left, bottom, width, height]
var_slider1 = Slider(
    ax=slider_ax1,
    label='Border 1',
    valmin=0,
    valmax=size_exg-10,
    valstep=1,
    valinit=initial_var1,
    track_color="white",
    color="0.5"
)

slider_ax2 = plt.axes([0.13, 0.1, 0.34, 0.03])  # [left, bottom, width, height]
var_slider2 = Slider(
    ax=slider_ax2,
    label='Border 2',
    valmin=0,
    valmax=size_exg-10,
    valstep=1,
    valinit=initial_var2,
    track_color="white",
    color="0.5"
)

rax = plt.axes([0.825, 0.1, 0.075, 0.075])
tax = plt.axes([0.55, 0.1, 0.075, 0.075])
radio_button = RadioButtons(rax,('linear', 'log'),activecolor="white")
transpose_button = Button(tax,"Transpose")

# The update function to call when the slider's value changes
def update(val):

    current_var1 = var_slider1.val  # Get the current slider value
    current_var2 = var_slider2.val

    diag1.set_xdata([-current_var1,size_exg-current_var1])
    diag2.set_xdata([-current_var2,size_exg- current_var2])
    tdat1.set_text(str(int(np.round(primary_Ex[max(current_var1,current_var2)])))+" keV")
    tdat2.set_text(str(int(np.round(primary_Ex[min(current_var1, current_var2)]))) + " keV")

    diag_arr = f(current_var1, current_var2)

    line.set_ydata(diag_arr)  # Recalculate and set new y-data
    line.set_xdata(primary_Ex[0:len(diag_arr)])
    ax2.set_ylim(5,np.max(diag_arr)*1.1)
    ax2.set_xlim(0,np.round(primary_Ex[size_exg-np.max([current_var1, current_var2])-1]))
    d_area.set_xy([(0,current_var2),(0,current_var1),(size_exg,size_exg+current_var1),(size_exg,size_exg+current_var2)])
    fig.canvas.draw_idle()  # Request a redraw of the figure

def clickfunc(label):
    current_var1 = var_slider1.val
    current_var2 = var_slider2.val

    ax2.set_yscale(label)
    ax2.set_ylim(5,np.max(f(current_var1, current_var2))*1.1)
    fig.canvas.draw_idle()

def bclickfunc(label):
    global primary_vals_rs
    current_var1 = var_slider1.val  # Get the current slider value
    current_var2 = var_slider2.val
    tmp=primary_vals_rs.transpose()
    primary_vals_rs=tmp
    primary_imshow.set_data(np.log(tmp+1))
    diag_arr = f(current_var1, current_var2)
    line.set_ydata(diag_arr)  # Recalculate and set new y-data
    line.set_xdata(primary_Ex[0:len(diag_arr)])
    ax2.set_ylim(5, np.max(diag_arr) * 1.1)
    ax2.set_xlim(0, np.round(primary_Ex[size_exg - np.max([current_var1, current_var2]) - 1]))
    d_area.set_xy([(0, current_var2), (0, current_var1), (size_exg, size_exg + current_var1),
                   (size_exg, size_exg + current_var2)])
    fig.canvas.draw_idle()
# Connect the update function to the slider

var_slider1.on_changed(update)
var_slider2.on_changed(update)
radio_button.on_clicked(clickfunc)
transpose_button.on_clicked(bclickfunc)

plt.show()
input()