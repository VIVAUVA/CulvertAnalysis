__author__ = 'Leo Cao'
from Tkinter import *
import tkFileDialog
from ProcessShapeFile import *


def askdirectory_culvert():
    fileName = tkFileDialog.askopenfilename()
    if fileName:
        str_var_culvert.set(fileName)


def askdirectory_shape():
    fileName = tkFileDialog.askopenfilename()
    if fileName:
        str_var_shape.set(fileName)


def askdirectory_output():
    fileName = tkFileDialog.asksaveasfilename()
    if fileName:
        str_var_output.set(fileName)


def UserFileInput(status, name, frame):
    optionFrame = frame
    optionLabel = Label(optionFrame, text=name)
    optionLabel.pack(side=LEFT)
    text = status
    str_var = StringVar(root)
    str_var.set(text)
    w = Entry(optionFrame, textvariable=str_var)
    w.pack(side=LEFT)
    optionFrame.pack()
    return w, str_var


def execute():

    print str_var_culvert.get(), str_var_shape.get(), str_var_output.get(), str_var_tri.get(), entry_radius.get()

    results_gui = extract_points_in_culvert(str_var_culvert.get(), str_var_shape.get(), str_var_output.get(),
                              tri_shape_file=str_var_tri.get(), input_search_radius=float(entry_radius.get()))

    culvert_info_display(results_gui)

if __name__ == '__main__':
    root = Tk()

    frame1 = Frame(root)
    w, str_var_culvert = UserFileInput("", "Culvert file", frame1)
    culvert_btn = Button(root, text='Culvert file (Excel)', command=askdirectory_culvert)
    culvert_btn.pack(side=TOP)
    frame1.pack(side=TOP)

    frame2 = Frame(root)
    w2, str_var_shape = UserFileInput("", "Shape File", frame2)
    shape_btn = Button(root, text='Shape file', command=askdirectory_shape)
    shape_btn.pack(side=TOP)
    frame2.pack(side=TOP)

    frame3 = Frame(root)
    w3, str_var_output = UserFileInput("", "Output File", frame3)
    frame3.pack()

    frame4 = Frame(root)
    w4, str_var_tri = UserFileInput("", "Triangulation output File", frame4)
    frame4.pack()

    L1 = Label(root, text="Search Radius for Culvert in feet (default value: 50)")
    entry_radius = Entry(root)
    L1.pack()
    entry_radius.pack()

    exe_btn = Button(root, text='Execute', command=execute)
    exe_btn.pack(side=BOTTOM)

    root.mainloop()