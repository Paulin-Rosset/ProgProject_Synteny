# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 15:35:05 2023

@author: user
"""
#GCF_001566615.1
#'WP_000241659.1'
#ex thread
#https://copyprogramming.com/howto/python-execute-python-script-threading-target-code-example

#python essais2.py


import tkinter as tk
import random as rd
import os
import sys
from threading import Thread
import module
import tkinter as tk
from PIL import Image, ImageTk
 
class StdoutRedirector(object):
    def __init__(self, text_widget):
        self.text_widget = text_widget
 
    def write(self, s):
        self.text_widget.insert('end', s)
        self.text_widget.see('end')

class Appli(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)
        self.height = 800
        self.width = 1400
        self.genome_var=tk.StringVar()
        self.protein_var=tk.StringVar()
        self.os_var=tk.StringVar()
        self.loading=False
        self.creer_widgets()

    def creer_widgets(self):
        # création canevas
        self.canv = tk.Canvas(self, bg="light gray", height=self.height,
                              width=self.width)
        self.canv.pack(side=tk.LEFT)
        # boutons
        self.genome = tk.Label(self, text="genome",font=('calibre',40, 'bold'))
        self.protein= tk.Label(self, text="protein",font=('calibre',40, 'bold'))
        self.os = tk.Label(self, text="OS",font=('calibre',40, 'bold'))
        self.genome_entry = tk.Entry(self,textvariable = self.genome_var, bd=5,font=('calibre',35,'normal'))
        self.protein_entry = tk.Entry(self,textvariable = self.protein_var, bd=5,font=('calibre',35,'normal'))
        self.os_entry = tk.Entry(self,textvariable = self.os_var, bd=5,font=('calibre',35,'normal'))
        self.submit_btn=tk.Button(self,text = 'Submit', height= 3, width=20,command = self.submit)
  
        self.genome.pack(side=tk.TOP)
        self.genome_entry.pack(side=tk.TOP)
        self.protein.pack(side=tk.TOP)
        self.protein_entry.pack(side=tk.TOP)
        self.os.pack(side=tk.TOP)
        self.os_entry.pack(side=tk.TOP)
        self.submit_btn.pack(side=tk.TOP)

        self.bouton_quitter = tk.Button(self, text="Quitter",
                                        command=self.quit)
        self.bouton_quitter.pack(side=tk.BOTTOM)
        
    def submit(self):
        print("executé")
        self.genome_use = self.genome_var.get()
        self.protein_use = self.protein_var.get()
        self.os_use = self.os_var.get()

        print("Genome chosen : " + self.genome_use)
        print("Protein chosen : " + self.protein_use)

        self.genome_var.set("")
        self.protein_var.set("")

        self.textbox = tk.Text(self.canv, wrap='word',height=45) #char pour couper au carac, word pour couper au mot
        self.textbox.pack()
        self.textbox.configure(font=("Times New Roman", 20, "italic"))
 
        ## Voici le truc, on pointe le sys.stdout vers le textbox
        sys.stdout = StdoutRedirector(self.textbox)

        thread = Thread(target = module.main(self.genome_use,self.protein_use,self.os_use))
        thread.start() # This code will execute in parallel to the current code
        thread.join()
 
        ## Ouverture du fichier
        image = Image.open('synt.png')
        image.show()


if __name__ == "__main__":
    app = Appli()
    app.title("BLAST Interface")
    app.mainloop()

