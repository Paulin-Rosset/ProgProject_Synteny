#Modules nécessaires:

import sys
import os
import pandas as pd
from Bio import SeqIO
import subprocess
from Bio.Blast import NCBIXML
import platform
from pygenomeviz import GenomeViz #Nécessité d'installer pygenomeviz: pip install pygenomeviz
from matplotlib.lines import Line2D
from datetime import datetime
from matplotlib.patches import Patch

##Fonctions de notre module:

#Fonctions à adapter pour run le code
def get_list_os():
    list_os = ["linux","win","macosx"]
    return list_os


#Fonctions partie 1 :
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
    dataframe.index = dataframe.index + 1
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
    return list_prot_info.head(1)

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

#Fonction partie 2:

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

#Fonctions partie 3 :

#Fonction lançant un blastp pour la protéine sur le génome indiqué:
def blastp(genome,record,my_os_type):
    """
    Description :
    ------------
        Fonction lançant un blastp pour la protéine sur le génome indiqué et créant un fichier xml à partir des résultats.
        
    Arguments :
    -----------
        genome : string , le génome d'intérêt.
        record : objet seq_record, informations sur la protéine étudiée.

    Returns :
    ---------
        Ne retourne rien, crée temporairement un fichier xml jusqu'à délétion dans une autre fonction.
    """
    #On crée le fichier fasta.
    SeqIO.write(record, "proteine.fasta", "fasta")
    
    #Sélection de l'OS.
    l_cmdos = ["data/Blast+/linux/bin/blastp","data/Blast+/win/bin/blastp.exe", "data/Blast+/macosx/bin/blastp"]
    cmdos = l_cmdos[get_list_os().index(my_os_type)]

    #On définit la commande à run par python.
    cmd_line = [cmdos,"-query","./proteine.fasta","-db","data/genomes/"+genome+"/"+genome,"-evalue","0.001","-out","result.xml","-outfmt", "5"]
    subprocess.run(cmd_line)

#Fonction récupérant le meilleur hit d'un blastp à partir d'un fichier xml contenant ses résultats.
def best_hit(file):
    """
    Description :
    ------------
        Fonction récupérant le meilleur hit, s'il existe, pour des résulats de blastp dans un tableau xml.
        
    Arguments :
    -----------
        file : string, le fichier xml contenant le résultat du blastp.

    Returns :
    ---------
        Retourne le meilleur hit sous forme d'un objet Alignement
    """
    
    with open(file) as handle:
        record = NCBIXML.read(handle)

    #Compteur.
    t = 0
    
    #Au cas où le gène n'est pas sur le gènome.
    result_hit = "Pas de résultat"
    
    for i in record.alignments:
        #Si t = 1 alors on a déja pris le meilleur hit qui est le premier.
        if t == 1:
            break
        else:
            t+=1
        result_hit = i

    #On supprimer le fichier xml du blast.
    os.remove("result.xml")

    #On supprime le fichier fasta crée plus tôt.
    os.remove("proteine.fasta")

    #On renvoie le résultat qui est un objet blast.
    return result_hit

#Fonctions partie 4 :

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
        Retourne une liste avec en premier l'ID de la protéine en amont et en second la protéine en aval.
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
    if recup_info_prot(genome,proteine)["Strand"][position_prot] == "+":
        proteine_aval = prot1
        proteine_amont = prot2
    else:
        proteine_aval = prot2
        proteine_amont = prot1

    return [proteine_amont, proteine_aval]

#Fonctions partie 6:

#Fonction associant -1 aux protéines sur le brin - et inversement.
def strand_protein(genome,dict):
    """
    Description :
    ------------
       Fonction associant -1 aux protéines sur le brin - et inversement.
        
    Arguments :
    -----------
        genome : string , le génome d'intérêt.
        dict : dictionary , dictionnaire regroupant pour chaque génome la protéine blastée associée.

    Returns :
    ---------
        Retourne 1 si protéine sur le brin + et -1 si protéine sur le brin -.
    """
    if dict[genome].Strand.item() == "-":
        return -1
    else : 
        return 1

#Fonction affichant les 3 protéines avec leur sens sur le brin.
def info_affichage_genome(genome, dict_aval, dict_amont, dict_centre):
    #Nom du brin, on donne le nom du génome.
    name = genome

    #Liste des positions des 3 protéines pour définir la taille de la visualisation.
    list_position = []
    
    #On gère si jamais la protéine centrale, en amont ou en aval n'a pas été conservée sur le génome. Si elle n'est pas conservée on a une chaine de caractères comme value du dictionnaire.
    if type(dict_amont[genome]) != str:
        #Comme la protéine est présente on compte sa position dans la liste des positions
        list_position += [dict_amont[genome].End.item(),dict_amont[genome].Begin.item()]
        
    if type(dict_aval[genome]) != str:
        list_position += [dict_aval[genome].End.item(),dict_aval[genome].Begin.item()]
        
    if type(dict_centre[genome]) != str:
        list_position += [dict_centre[genome].End.item(),dict_centre[genome].Begin.item()]

    #Si la liste des positions est vide on a aucun des gènes sur le génome blasté.
    if list_position == []:
        return {"name" : name, "size" : 0, "cds_list" : ((0,0,0,(0,0),""),(0,0,0,(0,0),""),(0,0,0,(0,0),""))}

    #Si on arrive ici alors la liste n'est pas nulle.
    #La fin de la zone de visualisation est le maximum de la liste des positions (avec une marge de 500pb).
    fin=max(list_position) + 500

    #Le début de la zone de visualisation est le minimum de la liste des positions (avec une marge de 500pb).
    debut=min(list_position) - 500

    #Taille de la zone de visualisation.
    taille = fin - debut 

   
    #On définit la taille de la protéine sur la zone de visualisation pour les 3 protéines.
    
    if type(dict_amont[genome]) == str:
        a = (0,0,0,(0,0),"")
    else:
        begin_cds_amont  = taille - (fin - dict_amont[genome].Begin.item())
        end_cds_amont    = taille - (fin - dict_amont[genome].End.item())
        signe_amont = strand_protein(genome, dict_amont)
        true_position_amont = (dict_amont[genome].Begin.item(),dict_amont[genome].End.item())
        protein_id_amont = dict_amont[genome].Protein_Id.item()
        a = (begin_cds_amont, end_cds_amont, signe_amont, true_position_amont, protein_id_amont)
    
    if type(dict_aval[genome]) == str:
        b = (0,0,0,(0,0),"")
    else:
        begin_cds_aval   = taille - (fin - dict_aval[genome].Begin.item())
        end_cds_aval     = taille - (fin - dict_aval[genome].End.item())
        signe_aval  = strand_protein(genome, dict_aval)
        true_position_aval = (dict_aval[genome].Begin.item(),dict_aval[genome].End.item())
        protein_id_aval = dict_aval[genome].Protein_Id.item()
        b = (begin_cds_aval, end_cds_aval, signe_aval, true_position_aval, protein_id_aval)
    
    if type(dict_centre[genome]) == str:
        c = (0,0,0,(0,0),"")
    else:
        begin_cds_centre = taille - (fin - dict_centre[genome].Begin.item())
        end_cds_centre   = taille - (fin - dict_centre[genome].End.item())
        signe_centre = strand_protein(genome, dict_centre)
        true_position_centre = (dict_centre[genome].Begin.item(),dict_centre[genome].End.item())
        protein_id_centre = dict_centre[genome].Protein_Id.item()
        c = (begin_cds_centre, end_cds_centre, signe_centre, true_position_centre, protein_id_centre)
    
    cds_list = (a, b, c)
    
    #On retourne les informations de visualisation.
    return {"name" : name, "size" : taille, "cds_list" : cds_list}

#Mane de notre module:

def main(genome,proteine,my_os_type):
    #Pour calculer le temps que met le script à s'exécuter:
    from datetime import datetime
    starttime = datetime.now()

    #Partie 0
    #Gérer le nombre d'arguments par utilisation fonction main
    var = ['',genome,proteine,my_os_type] 

    #Vérification des arguments en entrée, l'os est-il bien orthographié.
    list_os = ["linux","win","macosx"]

    if len(var) != 4:
        sys.exit("Fournir le bon nombre d'arguments svp.")
    else:
        if my_os_type not in get_list_os():
            sys.exit("Fournir un os correct : linux, win ou macosx")
            
    #Vérification des arguments en entrée, l'os est-il correct.
    result_platform = {'linux' : 'Linux' , 'macosx' : 'Darwin' , 'win' : 'Windows'} 
    
    if platform.system() != result_platform[my_os_type]:
        sys.exit(f"Fournir votre os, qui est : {platform.system()}")

    #Partie 1
    ##Si le nombre d'argument est bien de 2 et que le génome existe, 
    #alors on récupère les protéines présentes dans celui-ci et on
    #récupère les informations sur la protéine demandée.
    if genome in recup_genome("data/Ecoli_genomes_refseq.xlsx"):
        liste_proteines = recup_prot(genome)
        list_proteine_info = recup_info_prot(genome,proteine)
    else:
        sys.exit("Attention, le génome demandé n'est pas dans les 30 génomes donnés.")
        
    #Partie 2
    #On récupère la séquence de notre protéine d'intérêt dans notre génome d'intérêt avec quelques informations.
    record = recup_seq(genome,proteine)
    
    #Partie 3
    #On récupère les 29 autres génomes de la liste :
    list_genome_autre = recup_genome("data/Ecoli_genomes_refseq.xlsx")
    list_genome_autre.remove(genome)

    #Pour ne garder que l'ID.
    chars = 'ref|'
    
    #On crée un dictionnaire avec comme première valeur le génome de référence.
    dict_proteine_info = {genome : recup_info_prot(genome,proteine)}
    
    for genom in list_genome_autre:
        #On blast sur l'un des 29 génomes.
        blastp(genom,record,my_os_type)

        #Meilleur résultat.
        hit = best_hit("result.xml")

        #On gère si jamais le blastp n'avait pas eu de hit.
        if hit != "Pas de résultat":
            list_proteine_info = recup_info_prot(genom,hit.hit_id.translate(str.maketrans('', '', chars)))
            dict_proteine_info[genom] = list_proteine_info
        else:
            dict_proteine_info[genom] = "Pas de résultat"

    #Partie 4
    #On définit la protéine en amont et en aval.
    proteine_amont = amont_aval(genome,proteine)[0]
    proteine_aval = amont_aval(genome,proteine)[1]

    #Partie 5
    #Pour la protéine en amont, même procédé qu'en partie 3:
    record_amont = recup_seq(genome,proteine_amont)
    
    dict_proteine_amont_info = {genome : recup_info_prot(genome,proteine_amont)}
    
    for genom in list_genome_autre:
        #On blast sur l'un des 29 génomes.
        blastp(genom,record_amont,my_os_type)

        #Meilleur résultat.
        hit = best_hit("result.xml")

        #On gère si jamais le blastp n'avait pas eu de hit.
        if hit != "Pas de résultat":
            list_proteine_amont_info = recup_info_prot(genom,hit.hit_id.translate(str.maketrans('', '', chars)))
            dict_proteine_amont_info[genom] = list_proteine_amont_info
        else:
            dict_proteine_amont_info[genom] = "Pas de résultat"
    
    
    #Pour la proteine en aval, même procédé qu'avant.
    record_aval = recup_seq(genome,proteine_aval)
    
    dict_proteine_aval_info = {genome : recup_info_prot(genome,proteine_aval)}
    
    for genom in list_genome_autre:
        #On blast sur l'un des 29 génomes.
        blastp(genom,record_aval,my_os_type)

        #Meilleur résultat.
        hit = best_hit("result.xml")

        #On gère si jamais le blastp n'avait pas eu de hit.
        if hit != "Pas de résultat":
            list_proteine_aval_info = recup_info_prot(genom,hit.hit_id.translate(str.maketrans('', '', chars)))
            dict_proteine_aval_info[genom] = list_proteine_aval_info
        else:
            dict_proteine_aval_info[genom] = "Pas de résultat"
            

    #Partie 6
    #On crée les infos de visualtisation pour chaque génome.
    genome_visual=[]
    for genome_id in dict_proteine_info.keys() :
        genome_visual.append(info_affichage_genome(genome_id,dict_proteine_aval_info,dict_proteine_amont_info,dict_proteine_info))

    #On fait 3 boucles pour afficher 3 images.
    
    #Couleur en fonction de la place du gène, pour observer l'ordre.
    color = ["turquoise","plum","coral"]

    #Echelle.
    gv = GenomeViz(tick_style="axis")
    liste_trop_grand=[]
    #On place chaque gène sur le génome.
    for elem in genome_visual:
        name, size, cds_list = elem["name"], elem["size"], elem["cds_list"]
        if elem['size'] >= 20000 :
            track = gv.add_feature_track(elem["name"], taille)
            track.add_feature(500, 1000, 1, label = ("La distance entre les gènes est trop grande pour l'affichage"), linewidth = 1, labelrotation = 0, labelvpos = "top", labelhpos = "center", labelha = "center", labelsize = 14)
            liste_trop_grand.append([name, size, cds_list])
        elif elem['cds_list'] == ((0,0,0,(0,0),""),(0,0,0,(0,0),""),(0,0,0,(0,0),"")):
            track = gv.add_feature_track(elem["name"],taille)
            track.add_feature(500, 1000, 1, label = ("Les trois gènes sont absents du génome"), linewidth = 1, labelrotation = 0, labelvpos = "top", labelhpos = "center", labelha = "center", labelsize = 14)
        else :
            
            taille = size
            track = gv.add_feature_track(name, size)
            for idx, cds in enumerate(cds_list, 1):
                start, end, strand, position, protein_id = cds
                if position != (0,0):
                    if name == genome:
                        track.add_feature(0, 0, 1, label = "* Génome de référence *", labelrotation = 0, labelvpos = "top", labelhpos = "center", labelha = "center",labelsize = 8)
                    track.add_exon_feature([(start,end)], strand,exon_labels = [protein_id], labelrotation = 0, labelha = "center", exon_label_kws = {"y": 0, "va": "center", "color": "black"},labelsize = 9)
                    track.add_feature(start, end, strand, label = str(position[0]) + "             " + str(position[1]), linewidth = 1, labelrotation = 0, labelvpos = "top", labelhpos = "center", labelha = "center", facecolor = color[idx-1],labelsize = 10)

    #On gère l'affichage des légendes.
    
    #On trace la figure.
    fig = gv.plotfig()

    #On gère l'affichage des légendes.
    handles = [
    Patch(color = "coral", label = "Centre"),
    Patch(color = "turquoise", label = "Amont"),
    Patch(color = "plum", label = "Aval")
    ]
    
    fig.legend(handles = handles, bbox_to_anchor = (1, 1))  

    gv2 = GenomeViz(tick_style="axis")
    if liste_trop_grand!=[]:
        for i in liste_trop_grand :
                name2, size2, cds_list2=i[0],i[1],i[2]
                taille2 = size2
                track2 = gv2.add_feature_track(name2+"_2eme", size2)
                for idx, cds in enumerate(cds_list2, 1):
                    start, end, strand, position, protein_id = cds
                    if position != (0,0):
                        if idx == 1:
                            track2.add_feature(taille2/2, taille2/2, 1, label = "* Grande distance entre gènes *", labelrotation = 0, labelvpos = "top", labelhpos = "center", labelha = "center",labelsize = 8)
                        track2.add_exon_feature([(start,end)], strand,exon_labels = [protein_id], labelrotation = 0, labelha = "center", exon_label_kws = {"y": 0, "va": "center", "color": "black"},labelsize = 9)
                        track2.add_feature(start, end, strand, label = str(position[0]) + "             " + str(position[1]), linewidth = 1, labelrotation = 0, labelvpos = "top", labelhpos = "center", labelha = "center", facecolor = color[idx-1],labelsize = 10)
        fig2 = gv2.plotfig()
        #On gère l'affichage des légendes.
        handles = [
        Patch(color = "coral", label = "Centre"),
        Patch(color = "turquoise", label = "Amont"),
        Patch(color = "plum", label = "Aval")]
    
        fig2.legend(handles = handles, bbox_to_anchor = (1, 1))  
    else :
        track = gv2.add_feature_track(" ", 5000)
        track.add_feature(2500,2500,1,label="* Il n'y a aucune distance entre les gènes qui soit trop grande pour l'affichage *",labelrotation = 0, labelvpos = "top", labelhpos = "center", labelha = "center",labelsize = 8)
        fig2=gv2.plotfig()

    k = 1
    gv3 = GenomeViz(tick_style="axis")
    for elem in genome_visual:
        if elem["name"]==genome :
            name0, size0, cds_list0 = f'{elem["name"]} - 1', elem["size"], elem["cds_list"]
        name3, size3, cds_list3 = elem["name"], elem["size"], elem["cds_list"]
        if elem['cds_list'] == ((0,0,0,(0,0),""),(0,0,0,(0,0),""),(0,0,0,(0,0),"")):
            track = gv.add_feature_track(elem["name"] + ' ',taille)
            track.add_feature(500, 1000, 1, label = ("Les trois gènes sont absents du génome"), linewidth = 1, labelrotation = 0, labelvpos = "top", labelhpos = "center", labelha = "center", labelsize = 14)
        elif elem["size"]<=20000 :
            
            track3 = gv3.add_feature_track(name0, size0)
            for idx, cds in enumerate(cds_list0, 1):
                start, end, strand, position, protein_id = cds
                if position != (0,0):
                    track3.add_feature(0, 0, 1, label="* Genome de réference *", linewidth=1, labelrotation=0, labelvpos="top", labelhpos="center", labelha="center")
                    track3.add_exon_feature([(start,end)], strand,exon_labels = [protein_id], labelrotation = 0, labelha = "center", exon_label_kws = {"y": 0, "va": "center", "color": "black"},labelsize = 9)
                    track3.add_feature(start, end, strand, label = str(position[0]) + "             " + str(position[1]), linewidth = 1, labelrotation = 0, labelvpos = "top", labelhpos = "center", labelha = "center", facecolor = color[idx-1],labelsize = 10)
            
            track3 = gv3.add_feature_track(name3+"_3", size3)
            for idx, cds in enumerate(cds_list3, 1):
                start, end, strand, position, protein_id = cds
                if position != (0,0):
                    track3.add_exon_feature([(start,end)], strand,exon_labels = [protein_id], labelrotation = 0, labelha = "center", exon_label_kws = {"y": 0, "va": "center", "color": "black"},labelsize = 9)
                    track3.add_feature(start, end, strand, label = str(position[0]) + "             " + str(position[1]), linewidth = 1, labelrotation = 0, labelvpos = "top", labelhpos = "center", labelha = "center", facecolor = color[idx-1],labelsize = 10)

            # Add links between "genome référence" and "genome ième"
            if (cds_list3[0][0], cds_list3[0][1])!=(0,0):
                gv3.add_link((name0, cds_list0[0][0], cds_list0[0][1]), (name3+"_3", cds_list3[0][0], cds_list3[0][1]), normal_color="paleturquoise", inverted_color="lime", curve=True)
            if (cds_list3[1][0], cds_list3[1][1])!=(0,0):
                gv3.add_link((name0, cds_list0[1][0], cds_list0[1][1]), (name3+"_3", cds_list3[1][0], cds_list3[1][1]), normal_color="thistle", inverted_color="lime", curve=True)
            if (cds_list3[2][0], cds_list3[2][1])!=(0,0):
                gv3.add_link((name0, cds_list0[2][0], cds_list0[2][1]), (name3+"_3", cds_list3[2][0], cds_list3[2][1]), normal_color="lightcoral", inverted_color="lime", curve=True)
        k+=1
        name0 = genome + f' - {k}'

    fig3=gv3.plotfig()

    #On save la visualisation de la synténie.
    fig.savefig("synt1.png")
    fig2.savefig("synt2.png")
    fig3.savefig("synt3.png")
    
    #Pour connaître le temps à qu'a mis le script à run:
    difference = datetime.now() - starttime

    #Affichage des résultats:
    print(f"--------------------Protéine intérêt informations---------------------")
    print(f"\nVous avez lancé une analyse de synténie sur la protéine", proteine, "sur le génome d'E. coli", genome)
    print(f"\nLes informations sur notre protéine d'intéret sont : \n", list_proteine_info)
    print(f"\nLa séquence de la protéine d'intéret est : ", record.seq)
    print(f"\n-------------------Protéine amont/aval informations---------------------")
    print(f"\nLa protéine en amont est", proteine_amont)
    print(f"\nLes informations sur la protéine en amont sont :\n", recup_info_prot(genome,proteine_amont))
    print(f"\nLa protéine en aval est", proteine_aval)
    print(f"\nLes informations sur la protéine en amont sont :\n", recup_info_prot(genome,proteine_aval))
    print(f"\nLes blasts se sont effectués en",difference)
