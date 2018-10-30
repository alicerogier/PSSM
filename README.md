# PSSM

Ce repository contient :

*Un fichier matricePSSM.py qui contient quatre fonctions :
-matricePSSM(alignement, norm=False):
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
 -afficherMatricePSSM(matrice_PSSM, nligne, ncol)
    """
    Prend en entrée une matrice, son nombre de lignes et de colonnes.
    Affiche la matrice en entier
    """
 -remplissage_Dictionnaire_AA ()
    """
    Fonction qui construit un dictionnaire en attribuant à chaque lettre d'un acide aminé (clé) sa position dans la matrice PSSM (valeur)
    """
  -infos_Alignement(alignement)
    """
    Prend en entrée un fichier d'alignement
    Retourne la longueure de l'alignement et le nombre de séquences sous forme de tableau
    """
  Et une fonction main, qui permet de tester rapidement le code.  
*Deux fichiers différents (.malign) à mettre en arguments pour tester les code.
  
