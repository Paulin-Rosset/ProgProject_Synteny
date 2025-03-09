# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 15:35:06 2023

@author: user
"""
import tkinter as tk
from tkinter import ttk
import random as rd
import os
import sys
from threading import Thread
import module
import tkinter as tk
from PIL import Image


class StdoutRedirector(object):
    """
    Classe ajoutée pour pouvoir rediriger la sortie standard dans un widget de la fenêtre.
    """
    def __init__(self, text_widget):
        """
        Description: 
        Initialisation de la classe. On récupère en argument la Zone de texte dans lequel on veut afficher nos résultats d'analyse
        
        Argument: 
        text_widget: tk.Text
        """
        self.text_widget = text_widget

    def write(self, s):
        """
        Description:
        Inscrit la chaîne de caractère s à la fin du widget en argument de la classe
        
        Argument: 
        s : string
        """
        self.text_widget.insert('end', s)
        self.text_widget.see('end')
        
    def flush(self):
        """
        Description:
        Fonction primordiale pour l'execution sur Windows uniquement.
        """
        pass

class Appli(tk.Tk):
    """
    Crée un environnement où l'on définit les attributs de notre objet: une fenêtre graphique.
    """
    def __init__(self):
        """
        Description:
        Initialisation de la classe. On définit les attributs initiaux de notre fenêtre graphique.
        """
        tk.Tk.__init__(self)
        #Dimensions de la fenêtre graphique
        self.height = 650
        self.width = 800
        #Variables stockant la chaîne de caractères en temps réelle sur les textbars respectives
        self.genome_var=tk.StringVar()
        self.protein_var=tk.StringVar()
        #La fenêtre s'affiche avec les widgets 
        self.creer_widgets()
        #Messages de bienvenue
        self.canv.create_text((self.width/2,self.height/10),text="Bienvenue dans notre interface",font=('arial','30','bold'),fill='orange')
        self.canv.create_text((self.width/2,self.height/5),text="d'analyse de synténie",font=('arial','30','bold'),fill='orange')
        self.canv.create_text((self.width/2,2.2*self.height/5),text="Veuillez renseignez soigneusement",font=('arial','20','bold'),fill='black')
        self.canv.create_text((self.width/2,2.5*self.height/5),text="les paramètres de l'analyse",font=('arial','20','bold'),fill='black')
        self.canv.create_text((3.5*self.width/5,9.5*self.height/10),text="par M.Scalabrino, J.Baldous, P.Rosset",font=('arial','15','bold'),fill='black')

    def creer_widgets(self):
        """
        Description: 
        Création et mise en place des différents widgets de la fenêtre graphique
        """
        #Définition du canvas (zone colorée à gauche)
        self.canv = tk.Canvas(self, bg="#CDFAF6", height=self.height,width=self.width)
        self.canv.pack(side=tk.LEFT)

        #Définition des widgets demandant l'entrée des paramètres
        self.intro = tk.Label(self, text="Choix des paramètres:",font=('arial',25, 'bold'))
        self.genome = tk.Label(self, text="Genome",font=('arial',20, 'bold'))
        self.protein= tk.Label(self, text="Protein",font=('arial',20, 'bold'))
        self.os = tk.Label(self, text="OS",font=('arial',20, 'bold'))
        self.genome_entry = tk.Entry(self,textvariable = self.genome_var, bd=5,font=('arial',25,'normal'))
        self.protein_entry = tk.Entry(self,textvariable = self.protein_var, bd=5,font=('arial',25,'normal'))
        self.os_set=ttk.Combobox(self,values=module.get_list_os(),font=('arial',25,'normal'))
        self.os_set.current(0)
        self.submit_btn=tk.Button(self,text = 'Submit',bg="#CDFAF6", font=('arial',25),height= 2, width=10,command = self.submit)

        #pack() permet de placer proprement les widgets en fonction des autres. On pack dans l'ordre d'affichage de haut en bas
        self.intro.pack(side=tk.TOP)
        self.genome.pack(side=tk.TOP)
        self.genome_entry.pack(side=tk.TOP)
        self.protein.pack(side=tk.TOP)
        self.protein_entry.pack(side=tk.TOP)
        self.os.pack(side=tk.TOP)
        self.os_set.pack(side=tk.TOP)
        self.submit_btn.pack(side=tk.TOP)
        
    def submit(self):
        """
        Description:
        Bouton qui lance l'analyse de synténie et affiche les résultats. On s'assure de la validité des paramètres.
        """
        #On stocke les variables soumises par l'utilisateur
        self.genome_use = self.genome_var.get()
        self.protein_use = self.protein_var.get()
        self.os_use = self.os_set.get()

        #On s'assure de la validité des variables: génome est dans notre data et que la protéine est bien présente dans ce génome
        if self.genome_use in os.listdir('./data/genomes') and self.protein_use in module.recup_prot(self.genome_use):
            
            print("Genome chosen : " + self.genome_use)
            print("Protein chosen : " + self.protein_use)
            print("OS chosen : " + self.os_use)

            #On créé une zone de texte où l'on affichera les résultats de notre analyse
            self.textbox = tk.Text(self.canv, wrap='word',height=45,bg="#CDFAF6") #char pour couper au carac, word pour couper au mot
            self.textbox.pack()
            self.textbox.configure(font=("arial", 20, "italic"))
 
            #On redirige la sortie standard du script vers la textbox
            sys.stdout = StdoutRedirector(self.textbox)

            #On effectue l'analyse de synténie en parallèle de la fenêtre graphique
            thread = Thread(target = module.main(self.genome_use,self.protein_use,self.os_use))
            thread.start() 


            #On ouvre le fichier illustrant notre analyse de synténie 
            image1 = Image.open('synt1.png')
            image2 = Image.open('synt2.png')
            image3 = Image.open('synt3.png')
            image1.show()
            image2.show()
            image3.show()
        else:
            self.canv.create_text((400,self.3.5*height/5),text="Au moins un des arguments est incorrect",font=('arial','25','bold'),fill='red')




if __name__ == "__main__":
    app = Appli()
    app.title("Synteny Analysis Interface")
    app.mainloop()
