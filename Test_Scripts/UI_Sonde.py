import tkinter as tk
from tkinter import ttk
#from tkinter.filedialog import askopenfilename
from tkinter import filedialog
from tkinter import *

import os


def import_file():
    global sounding_file_path
    file_path = filedialog.askopenfilename(title="Select a file", filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
    if file_path:
        #print("Selected file:", file_path)
        sounding_file_path = file_path


def submit_callback(e1, e2, sounding_file_path):
#     print("User entered : " + e1.get())
#     print("User entered : " + e2.get())
    #print("Will now run file with : " + e1.get() + e2.get() + sounding_file_path)
    os.system(f'python UI_Sonde_run.py --inputlat={e1.get()} --inputlon={e2.get()} --inputfile={sounding_file_path}')
    return None


# Creating tkinter window
window = tk.Tk()
window.title('UI-Sonde')
window.geometry('500x550')

input_label = tk.Label(window, text="Welcome to UI-Sonde!",
                       font = ("Ariel", 18), foreground="blue").grid(column = 0,row = 1, padx = 0, pady = 25)
description_label = tk.Label(window, text="Select a file and input station lat/lon to generate a Skew-T",
                             foreground="black", font = ("Ariel", 10)).grid(column = 0,row = 2, padx = 10, pady = 25)


# Input values
Label(window, text='Station Latitude (degrees)').grid(row=5)
Label(window, text='Station Longitude (degrees)').grid(row=6)
e1 = Entry(window)
e2 = Entry(window)
e1.grid(row=5, column=1)
e2.grid(row=6, column=1)

submit_button = tk.Button(window, text = "Submit" , command = lambda: submit_callback(e1, e2, sounding_file_path))
submit_button.grid(column=1, row=10)

# Create an "Import File" button
import_button = tk.Button(window, text="Import File", command=import_file)
import_button.grid(row=4, column=1)


mainloop()