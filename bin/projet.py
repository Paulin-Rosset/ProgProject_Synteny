import sys
import os
import pandas as pd
from Bio import SeqIO

##Notre dossier de travail est celui juste avant le dossier data.

#Récupération des données (protéine et génome)
var = sys.argv
genome = var[1]
proteine = var[2]
my_os_type = var[3]

#Vérification des arguments en entrée:
list_os = ["linux","win","macosx"]

if len(var) != 4:
    sys.exit("Fournir le bon nombre d'arguments svp.")
else:
    if my_os_type not in list_os :
        sys.exit("Fournir un os correct : linux, win ou macosx")


###PARTIE 1:

##Fonctions :
#Fonction créeant un dataframe à partir du fichier d'annotation du génome spécifié:
def recup_data(genome):
    """
    Description :
    ------------
        Fonction récupérant le dataframe du fichier d'annotation à du génome donnée en argument.

    Arguments :
    -----------
        genome : string

    Returns :
    ---------
        Retourne le dataframe correspondant.

    """
    ##On récupére en dataframe le fichier tsv.
    dataframe = pd.read_csv(f'data/genomes/{genome}/annotation_{genome}.tsv', sep = '\t')
    dataframe.index += 1
    return dataframe


#Fonction récupérant les informations sur la protéine du génome spécifié:
def recup_info_prot(genome,proteine):
    """
    Description :
    ------------
        Fonction récupérant les informations associée à la protéine et au génome rentrés en argument.

    Arguments :
    -----------
        genome : string , le génome d'intérêt.
        proteine : string , la protéine d'intérêt.

    Returns :
    ---------
        Retourne un dataframe correspondant à la ligne de la protéine du fichier d'annotation.
    """
    ##La gestion d'erreur sert si on exécute le code dans le mauvais dossier de travail ou que le dossier contenant le génome n'existe pas.
    try:
        ##On récupére en dataframe le fichier tsv.
        dataframe = recup_data(genome)
        ##Si la protéine qu'on a demandé n'est pas dans le génome on break le script, sinon on récupère en dictionaire les infos sur la protéine. A OPTIMISER
        if proteine in set(dataframe["Protein_Id"]):
            list_prot_info = dataframe[dataframe["Protein_Id"] == proteine]
        else:
            sys.exit("La protéine n'est pas contenu dans le génome donné.")
    except FileNotFoundError:
        print("Rentrez un argument correct svp.")
    return list_prot_info

#Fonction récupérant la liste des protéines dans le génome spécifié:
def recup_prot(genome):
    """
    Description :
    ------------
        Fonction récupérant la liste des protéines disponibles dans le génome d'intérêt.
        
    Arguments :
    -----------
        genome : string , le génome d'intérêt.
        
    Returns :
    ---------
        Retourne une liste avec tous les id des protéines du génome d'intérêt.
    """
    ##La gestion d'erreur sert si on exécute le code dans le mauvais dossier de travail ou que le dossier contenant le génome n'existe pas.
    try:
        ##On récupére en dataframe le fichier tsv.
        dataframe = recup_data(genome)
        ##On crée une liste avec toutes les protéines.
        list_prot = list(dataframe.Protein_Id)
    except FileNotFoundError:
        print("Rentrez un argument correct svp.")
    return list_prot


#Fonction récupérant une liste avec les 30 génomes:
def recup_genome(file):
    """
    Description :
    ------------
        Fonction récupérant la liste des 30 génomes par leur Assembly accesion.
        
    Arguments :
    -----------
        file : string , chemin d'accès au fichier, généralement "data/Ecoli_genomes_refseq.xlsx"

    Returns :
    ---------
        Retourne une liste composée des 30 génomes de notre projet.
    """
    ##La gestion d'erreur sert si on exécute le code dans le mauvais dossier de travail.
    try :
        ##On converit le fichier excell en argument en dataframe.
        classeur = pd.read_excel(file)
        ##On récupère seulement les génomes, représentés par le Assembly accesion.
        genom = list(classeur["Assembly Accession"])
    except FileNotFoundError:
        print("Rentrez un argument correct svp.")
    return genom

##Main:

##Si le nombre d'argument est bien de 2 et que le génome existe, 
##alors on récupère les protéines présentes dans celui-ci et on
##récupère les informations sur la protéine demandée.
if genome in recup_genome("data/Ecoli_genomes_refseq.xlsx"):
    liste_proteines = recup_prot(genome)
    list_proteine_info = recup_info_prot(genome,proteine)
else:
    sys.exit("Attention, le génome demandé n'est pas dans les 30 génomes donnés.")
    
#Vérification.
#print("La liste des différentes protéines du génome est : ", liste_proteines)
print("Les informations sur notre protéine d'intéret sont : ", list_proteine_info)


###PARTIE 2:
#Si on arrive à cette partie c'est que la protéine existe bien dans le génome et qu'on a récupéré les informations dans dict_proteine_info.

##Fonctions:
#Fonction qui récupère la séquence de la protéine du fichier fasta contenant les séquences de toutes les protéines du génome :
def recup_seq(genome,protein):
    """
    Description :
    ------------
        Fonction récupérant l'objet seq_record associé à notre protéine.
        
    Arguments :
    -----------
        genome : string , le génome d'intérêt.
        proteine : string , la protéine d'intérêt.

    Returns :
    ---------
        Retourne un objet seq_record contenant la séquence de notre protéine d'intérêt et son identifiant.

    """

    for seq_record in SeqIO.parse(f'data/genomes/{genome}/protein.faa',"fasta"):
        if seq_record.id == protein:
            return(seq_record)

##Main:
#On récupère la séquence de notre protéine d'intérêt dans notre génome d'intérêt.
sequence = recup_seq(genome,proteine).seq

#Vérification.
print("La séquence de la protéine d'intéret est : ", sequence)


###PARTIE 3: 

##Fonctions:


##Main:

#On récupère les 29 autres génomes de la liste :
list_genome_autre = recup_genome("data/Ecoli_genomes_refseq.xlsx")
list_genome_autre.remove(genome)


###PARTIE 4:

##Fonctions:
#Fonction cherchant la protéine en amont et en aval de la protéine donnée:

def amont_aval(genome,proteine):
    """
    Description :
    ------------
        Fonction cherchant la protéine en amont et en aval de notre protéine d'intérêt dans notre génome d'intérêt.
        
    Arguments :
    -----------
        genome : string , le génome d'intérêt.
        proteine : string , la protéine d'intérêt.

    Returns :
    ---------
        Retourne une liste avec en premier la protéine en amont et en second la protéine en aval.
    """
    
    #On récupére en dataframe le fichier tsv.
    dataframe = recup_data(genome)

    #On récupère la position dans le dataframe de notre protéine.
    position_prot = dataframe[dataframe["Protein_Id"] == proteine].index[0]

    #On fait une boucle pour gérer si la protéine est au début du génome ou à la fin (dans le tableau tsv).
    if position_prot == 1:
        print("Attention votre protéine est sur le premier gène du génome, comme il est circulaire ce n'est pas pénalisant")
        prot1 = dataframe.loc[position_prot + 1]["Protein_Id"]
        prot2 = dataframe.loc[dataframe.shape[0]]["Protein_Id"]
    elif position_prot == dataframe.shape[0]:
        print("Attention votre protéine est sur le dernier gène du génome, comme il est circulaire ce n'est pas pénalisant")
        prot1 = dataframe.loc[1]["Protein_Id"]
        prot2 = dataframe.loc[position_prot - 1]["Protein_Id"]
    else:
        prot1 = dataframe.loc[position_prot + 1]["Protein_Id"]
        prot2 = dataframe.loc[position_prot - 1]["Protein_Id"]

    #On gère le fait qu'un gène en brin + ou - n'a pas les mêmes définitions de aval et amont.
    if list_proteine_info["Strand"][position_prot] == "+":
        proteine_aval = prot1
        proteine_amont = prot2
    else:
        proteine_aval = prot2
        proteine_amont = prot1

    return [proteine_amont, proteine_aval]

##Main:

#On définit la protéine en amont et en aval.
proteine_amont = amont_aval(genome,proteine)[0]
proteine_aval = amont_aval(genome,proteine)[1]

#Vérification.
print("La protéine en amont est", proteine_amont)
print("La protéine en aval est", proteine_aval)


###PARTIE 5









