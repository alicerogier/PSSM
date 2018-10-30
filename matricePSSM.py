#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 10:06:04 2018

@author: alice
"""

import numpy as np

def matricePSSM(alignement, norm=False):
    """
    Prend en entrée un fichier d'alignement MUSCLE (.align).
    Retourne une matrice PSSM .
    La matrice PSSM est une matrice de dimension 20*longueure de l'alignement.
    Elle est de la forme : 
        
    position :          0   1   2   3   4   5  ... len(alignement)
    
    ----------------
    ordre des
    aa dans 
    matrice :

    0             A
    1             C
    2             D
    3             E
    4             F
    5             G
    6             H
    7             I
    8             K
    9             L
    10            M
    11            N
    12            P
    13            Q
    14            R
    15            S
    16            T
    17            V
    18            W
    19            Y
        
    """
    
    import sys 
    # stores PDB
    try :
        f = open(alignement, "r")
        
    except :
        print("L'ouverture a echoue")
        sys.exit(1)
        
    #lecture du fichier :
    
    lines = f.readlines()
    f.close()
    
    #variables:
    
    dico_AA_matrice=remplissage_Dictionnaire_AA()  #dictionnaire qui prend en clef la lettre de l'acide aminé, et en valeur la ligne correspondante dans la matrice PSSM
    nb_seq=infos_Alignement(alignement)[0] #nombre de séquences du fichier afin de normaliser la matrice si norm=True
    position=0 #position dans la séquence et dans la colonne de la matrice PSSM
    longueure_alignement=infos_Alignement(alignement)[1]#longueur de l'alignement
    matrice_PSSM=np.zeros((20,longueure_alignement)) #matrice de 0 de longueure 20*la longueure de l'alignement
    
    #Remplissage de la matrice PSSM :
    
    for line in lines :
        if line[0]==">" : #Si on est sur une ligne où il y a le nom de la séquence
            position=0 #On remet à 0 le compteur de position dans la séquence
        else : #Si on est sur une ligne de séquence
            for i in line : #On est au sein d'une séquence, i prend a pour valeur un char (lettre d'acide aminé)
                if i != "\n" and i != "-": #Si i ne correpond pas à un retour à la ligne, ni à un gap (il corrpond donc forcément à un acide aminé
                    matrice_PSSM[dico_AA_matrice[i],position]+=1 #Pour l'acide aminé considéré, à la position considérée, on ajoute 1
                    position+=1 #On incrémente le compteur de position au sein de la séquence considérée                                   
    if norm :#Si l'option norm=True on obtient des fréquences d'acides aminés à la place d'entiers
        for i in range (20):
            for j in range (longueure_alignement):
                matrice_PSSM[i][j]=matrice_PSSM[i][j]/nb_seq
        
    return(matrice_PSSM)
 
        
def afficherMatricePSSM(matrice_PSSM, nligne, ncol):
    """
    Prend en entrée une matrice, son nombre de lignes et de colonnes.
    Affiche la matrice en entier
    """
    for i in range (nligne):
        for j in range (ncol):
            print(matrice_PSSM[i][j]," ", end='')
        print("\n")              
                
def remplissage_Dictionnaire_AA () :
    """
    Fonction qui construit un dictionnaire en attribuant à chaque lettre d'un acide aminé (clé) sa position dans la matrice PSSM (valeur)
    """
    dico_AA_matrice={} #dictionnaire qui prend en clef la lettre de l'acide aminé, et en valeur la ligne correspondante dans la matrice PSSM
    liste_AA=["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"] #liste des acides aminés dans l'ordre des lignes de la matrice PSSM
    
    #attribution des couples clé-valeur du dicrtionnaire :
    k=0;
    
    for j in liste_AA:
        dico_AA_matrice[j]=k
        k+=1;
    
    return(dico_AA_matrice)
    
def infos_Alignement(alignement):
    """
    Prend en entrée un fichier d'alignement
    Retourne la longueure de l'alignement et le nombre de séquences sous forme de tableau
    """
    import sys 
    # stores PDB
    try :
        f = open(alignement, "r")
        
    except :
        print("L'ouverture a echoue")
        sys.exit(1)

    #lecture du fichier :
    
    lines = f.readlines()
    f.close()
    
    #variables :
    nb_seq=0
    tab_info=[0,0]
    seq=""
    
    for line in lines :
        if line[0]==">":
            #print(seq)
            #print(len(seq))
            nb_seq+=1
            seq=""
        else :
            for i in line:
                seq=seq+i
                
            
    tab_info[0]=nb_seq
    tab_info[1]=len(seq)
    
    return(tab_info)
    
    
def main() :
    from matricePSSM import matricePSSM
    from matricePSSM import infos_Alignement
    from matricePSSM import remplissage_Dictionnaire_AA
    from matricePSSM import afficherMatricePSSM
    matrice=matricePSSM("testx.malign",norm=True)
    afficherMatricePSSM(matrice,20,infos_Alignement("testx.malign")[1])
    
    
    
    
    
    
    
               
        
    
    