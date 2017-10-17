#----------------------------------------------------#
#--Author: Harrys Kon (Charalambos Konstantinou)-----#
#--W: https://harrys.fyi/----------------------------#
#--E: konharrys@gmail.com----------------------------#
#----------------------------------------------------#
#!/usr/bin/env python
import fileinput
import shutil
from tkinter import *
import time
root=Tk()



variable1=StringVar()
variable2=StringVar()
variable3=StringVar()
variable4=StringVar()
variable5=StringVar()
variable6=StringVar()
variable7=StringVar()
variable8=StringVar()
variable9=StringVar()
variable10=StringVar()
variable11=StringVar()
variable12=StringVar()
variable13=StringVar()
variable14=StringVar()
    
        
def update_label(filename):
    f = open(filename)
    lines = f.readlines()
    f.close()
    i = 0
    j= 0
    k = 0
    N=14
    while i<(len(lines)/N):
        state = []
        j = 0
        while (j<N):
            state.append(lines[k].strip().split())
            j = j + 1
            k = k + 1
        i = i + 1
        variable1.set("V" + str(state[0][0]) + " = " + str(state[0][1]) + "\t" + "\u03B8" + str(state[0][0]) + " = " + str(state[0][2]))
        variable2.set("V" + str(state[1][0]) + " = " + str(state[1][1]) + "\t" + "\u03B8" + str(state[1][0]) + " = " + str(state[1][2]))
        variable3.set("V" + str(state[2][0]) + " = " + str(state[2][1]) + "\t" + "\u03B8" + str(state[2][0]) + " = " + str(state[2][2]))
        variable4.set("V" + str(state[3][0]) + " = " + str(state[3][1]) + "\t" + "\u03B8" + str(state[3][0]) + " = " + str(state[3][2]))
        variable5.set("V" + str(state[4][0]) + " = " + str(state[4][1]) + "\t" + "\u03B8" + str(state[4][0]) + " = " + str(state[4][2]))
        variable6.set("V" + str(state[5][0]) + " = " + str(state[5][1]) + "\t" + "\u03B8" + str(state[5][0]) + " = " + str(state[5][2]))
        variable7.set("V" + str(state[6][0]) + " = " + str(state[6][1]) + "\t" + "\u03B8" + str(state[6][0]) + " = " + str(state[6][2]))
        variable8.set("V" + str(state[7][0]) + " = " + str(state[7][1]) + "\t" + "\u03B8" + str(state[7][0]) + " = " + str(state[7][2]))
        variable9.set("V" + str(state[8][0]) + " = " + str(state[8][1]) + "\t" + "\u03B8" + str(state[8][0]) + " = " + str(state[8][2]))
        variable10.set("V" + str(state[9][0]) + " = " + str(state[9][1]) + "\t" + "\u03B8" + str(state[9][0]) + " = " + str(state[9][2]))
        variable11.set("V" + str(state[10][0]) + " = " + str(state[10][1]) + "\t" + "\u03B8" + str(state[10][0]) + " = " + str(state[10][2]))
        variable12.set("V" + str(state[11][0]) + " = " + str(state[11][1]) + "\t" + "\u03B8" + str(state[11][0]) + " = " + str(state[11][2]))
        variable13.set("V" + str(state[12][0]) + " = " + str(state[12][1]) + "\t" + "\u03B8" + str(state[12][0]) + " = " + str(state[12][2]))        
        variable14.set("V" + str(state[13][0]) + " = " + str(state[13][1]) + "\t" + "\u03B8" + str(state[13][0]) + " = " + str(state[13][2]))
        root.update()
        time.sleep(0.5)



root.title("State Estimation Results")
label1=Label(root,fg="dark green",textvariable=variable1)
label1.pack()
label2=Label(root,fg="dark green",textvariable=variable2)
label2.pack()
label3=Label(root,fg="dark green",textvariable=variable3)
label3.pack()
label4=Label(root,fg="dark green",textvariable=variable4)
label4.pack()
label5=Label(root,fg="dark green",textvariable=variable5)
label5.pack()
label6=Label(root,fg="dark green",textvariable=variable6)
label6.pack()
label7=Label(root,fg="dark green",textvariable=variable7)
label7.pack()
label8=Label(root,fg="dark green",textvariable=variable8)
label8.pack()
label9=Label(root,fg="dark green",textvariable=variable9)
label9.pack()
label10=Label(root,fg="dark green",textvariable=variable10)
label10.pack()
label11=Label(root,fg="dark green",textvariable=variable11)
label11.pack()
label12=Label(root,fg="dark green",textvariable=variable12)
label12.pack()
label13=Label(root,fg="dark green",textvariable=variable13)
label13.pack()
label14=Label(root,fg="dark green",textvariable=variable14)
label14.pack()
stop_button=Button(root,text="stop",width=25,command=root.destroy)
stop_button.pack()
update_label("../files/outputSE.txt")
root.mainloop()









