from math import *
import matplotlib.pyplot as plt
import matplotlib.animation
import numpy as np
import copy as cp

## Donnees en entree
xCible = 0 #coordonees de la cible donnees par le traitement d'image en m
yCible = 0

xBateau = 0 #coordonees du bateau donnees par le traitement d'image en m
yBateau=0

thetaBateau = 0 #angle en degre du bateau donne par sa boussole numerique, par rapport a une reference

deltaT = 1 #intervalle de temps en secondes pour lequel les valeurs de position s'actualisent (avec le traitement d'image)

## Constantes a choisir
rayonRalentissement = 00 #rayon du cercle en cm dans lequel on diminue la vitesse du bateau pour ajuster l'angle
                          #on appelle dans la suite le cercle ayant pour centre la cible et pour rayon rayonRalentissement la "zone de ralentissement"
rayonArrivee = 10
puissanceHeliceConstante = 100 #la puissance qu'on "propose" au bateau lorsqu'il n'est pas dans le cercle min, entre 0% et 100%
puissanceHeliceMin = 20 #la puissance minimum pour l'hélice du bateau telle que le bateau arrive à entrer en mouvement
puissanceHeliceGrandAngle = 50 #la puissance qu'on "propose" au bateau lorsqu'il n'est pas dans le cercle min et qu'il doit effectuer un grand virage, entre 0% et 100%
grandAngle = 20 #definition d'un "grand angle" en degres

vitesseCible = 1 #vitesse cible du bateau lorsqu'il n'est pas dans le cercle min, en m/s


Kp = 1 #parametres du PID en angle
Ki = 0
Kd = 1

Lp = 1#parametres du PID en vitesse
Li = 1
Ld = 1


## PID et Asservissement
sommeErreursAngle = 0
erreurPrecedenteAngle = 0.001

sommeErreursVitesse = 0
erreurPrecedenteVitesse = 0.001

distancePrecedente = 0 #idéalement, il faudrait initialiser avec la distance initiale du bateau à (0,0)

# Remarque : on considere avoir fixe une origine des angles theta0. On notera par la suite :
#                   - thetaBateau l'angle du bateau par rapport a theta0
#                   - thetaCible l'angle fait par la droite passant par le bateau et la cible par rapport a theta0
#                   - thetaCommande l'angle que l'on ajoute à thetaBateau pour imposer une consigne au gouvernail du bateau par rapport a theta0
# Idee de l'asservissement :
#   On va asservir uniquement la direction du bateau a l'aide du PID
#   La vitesse du bateau, qu'on peut assimiler a la puissance de l'helice est constante au cours de la trajectoire, sauf quand le bateau entre dans la zone de ralentissement
#   Dans la zone de ralentissement, on ralentit la puissance de l'helice lineairement, donc la vitesse, et on continue d'asservir l'angle

def PIDangle (thetaCible , thetaBateau ): #retourne thetaCommande necessaire pour l'asservissement
    global sommeErreursAngle
    global erreurPrecedenteAngle
    if (thetaCible - thetaBateau)%360 <= 180:
        erreur = (thetaCible - thetaBateau)%360
    else:
        erreur = (thetaCible - thetaBateau)%360 -360
    sommeErreursAngle += erreur #integration
    sommeErreursAngle = sommeErreursAngle
    variationErreur = (erreur - erreurPrecedenteAngle) #derivation
    thetaCommande = (Kp * erreur + Ki * sommeErreursAngle + Kd * variationErreur) #formule du PID
    erreurPrecedenteAngle = erreur

    # if (thetaBateau-thetaCommande)%360 > 80 and (thetaBateau-thetaCommande)%360 < 280:
    #     if (thetaBateau-thetaCommande)%360 < 180:                    # Ce bloc if else permet d'éviter d'imposer en commande plus d'un demi tour,
    #         thetaCommande = (thetaBateau -80)                 # ce qui n'aurait pas de sens car le bateau tournerait dans le sens inverse a celui voulu
    #     else :
    #         thetaCommande = (thetaBateau +80)
    if thetaCommande>0:
        thetaCommande  = min(thetaCommande, 60)
    if thetaCommande<=0:
        thetaCommande  = max(thetaCommande, -60)

    return thetaCommande

def PIDvitesse (vitesseCible, vitesseBateau):#retourne vitesseCommande necessaire pour l'asservissement
    global sommeErreursVitesse
    global erreurPrecedenteVitesse
    erreur = vitesseCible - vitesseBateau
    sommeErreursVitesse += erreur #integration
    variationErreur = erreur - erreurPrecedenteVitesse #derivation
    vitesseCommande = Lp * erreur + Li * sommeErreursVitesse + Ld * variationErreur #formule du PID
    erreurPrecedenteVitesse = erreur
    return vitesseCommande


def vitesse (xBateau , yBateau): #vitesse du bateau en m/s
    global distancePrecedente
    distance = xBateau**2 + yBateau**2 #distance du bateau par rapport au point fixe arbitraire de coordonnees (0,0)
    deltaR = distance - distancePrecedente
    vitesseBateau = deltaR/deltaT
    distancePrecedente = distance
    if (distancePrecedente == 0) : #uniquement lorsque distancePrecedente n'est pas initialisee, a la premiere utilisation de la fonction
        return 0
    else :
        return vitesseBateau

def puissanceHelice (xCible , yCible , xBateau , yBateau , thetaBateau): #renvoie la puissance de l'hélice du bateau entre 0 et 100%
    distanceBateauCible, thetaCible = trajetCible (xCible , yCible , xBateau , yBateau , thetaBateau)
    vitesseBateau = vitesse(xBateau, yBateau)
    vitesseCommande = PIDvitesse (vitesseCible, vitesseBateau) + vitesseBateau
    #print('Vitesse Commande : '+str(vitesseCommande))

    if (distanceBateauCible > rayonRalentissement) : #si on est en dehors de la zone de ralentissement, on asservit la vitesse selon l'angle
        if (vitesseCommande > vitesseCible): #cas ou la vitesse est trop elevee
            return puissanceHeliceMin
        else : #cas ou la vitesse est trop faible
            if (abs(thetaBateau-thetaCible)>grandAngle):#cas ou l'angle est eleve
                return puissanceHeliceGrandAngle
            else :#cas ou l'angle est raisonnable
                return puissanceHeliceConstante
    else : #si on est dans la zone de ralentissement, on ralentit la vitesse lineairement et on continue d'asservir l'angle
        puissanceHeliceZone = max((distanceBateauCible/rayonRalentissement) * puissanceHeliceConstante, puissanceHeliceMin) #diminue la puissance de l'helice linéairement de façon a ajuster l'angle precisement
        puissanceHeliceZoneGrandAngle = max((distanceBateauCible/rayonRalentissement) * puissanceHeliceGrandAngle, puissanceHeliceMin)
        if (vitesseCommande > vitesseCible): #cas ou la vitesse est trop elevee
            return 0
        else : #cas ou la vitesse est trop faible
            if (abs(thetaBateau-thetaCible)>grandAngle):#cas ou l'angle est eleve
                return puissanceHeliceZoneGrandAngle
            else :#cas ou l'angle est raisonnable
                return puissanceHeliceZone



def trajetCible (xCible , yCible , xBateau , yBateau , thetaBateau):
    distanceBateauCible = sqrt((xCible-xBateau)**2 + (yCible-yBateau)**2) #distance entre le bateau et la cible
    # if (xCible == xBateau): #permet d'eviter la division par 0 dans l'arctangente
    #     if (yBateau>=yCible):
    #         thetaCible = 270
    #     else :
    #         thetaCible = 90
    # else :
    thetaCible = ((180/pi)*atan2((yCible-yBateau),(xCible-xBateau)))%360
    return distanceBateauCible, thetaCible


def asservissement(xCible , yCible , xBateau , yBateau , thetaBateau): #retourne l'angle (absolu) de commande des gouvernails et le pourcentage de puissance de l'helice (0 à 100)
    distanceBateauCible, thetaCible = trajetCible (xCible , yCible , xBateau , yBateau , thetaBateau)
    thetaCommande = PIDangle (thetaCible, thetaBateau)
    #thetaAbsolu = (thetaBateau + thetaCommande)%360
    puissanceHeliceAbsolue = puissanceHelice(xCible , yCible , xBateau , yBateau , thetaBateau)
    return thetaCommande, puissanceHeliceAbsolue

##Navigation : choix de l'ordre des cibles
limite = 0.7

def courbure (xb,yb,xc,yc,theta): #Cette fonction renvoie la courbure de la trajectoire qui est un arc de cercle entre le bateau et la cible. La courbure est alors egale a l'inverse du rayon du cercle

    if theta!=0:

        x0=(((yb-yc)*(yb+(xb/tan(theta))))+((xc**2-xb**2)/2)+((yc**2-yb**2)/2))/(xc-xb-((yc-yb)/tan(theta)))
        y0=(-x0/tan(theta))+(yb+(xb/tan(theta))) #x0 et y0 sont les coordonnees du centre du cercle
        R=sqrt((xb-x0)**2+(yb-y0)**2) #R est le rayon du cercle


    else:
        theta = 0.0001
        x0=(((yb-yc)*(yb+(xb/tan(theta))))+((xc**2-xb**2)/2)+((yc**2-yb**2)/2))/(xc-xb-((yc-yb)/tan(theta)))
        y0=(-x0/tan(theta))+(yb+(xb/tan(theta))) #x0 et y0 sont les coordonnees du centre du cercle
        R=sqrt((xb-x0)**2+(yb-y0)**2)


    return 1/R


def listeCiblesAcessibles (listeDesCibles,xb,yb,theta): #Fait la liste des cibles accessibles par le bateau
    listeCiblesAccessibles = []
    for i in range (len(listeDesCibles)-1) :
        if (courbure(xb,yb,listeDesCibles[i][0],listeDesCibles[i][1],theta)<= limite): #En considerant que la liste des cibles est un tableau de taille (le nombre de cibles)x2
            listeCiblesAccessibles.append(listeDesCibles[i])
    return listeCiblesAccessibles


def calculDistance (xb,yb,xc,yc): #Fonction qui calcule la distance entre la cible et le bateau

    distance = sqrt ((xb-xc)**2+(yb-yc)**2)
    return distance

def triCiblesDistance (xb,yb,listeCiblesAccessibles): #Cette fonction trie la liste des cibles accessibles par le bateau par ordre croissant de leur distance
    listeDistances = [] #tableau comportant la distance entre le bateau et chaque cible ainsi que le numero de la cible consideree
    listeCiblesTriees=[]
    for i in range (len(listeCiblesAccessibles)):
        distance = calculDistance (xb,yb,listeCiblesAccessibles[i][0],listeCiblesAccessibles[i][1])
        listeDistances.append((distance,i))
    listeDistances.sort()
    for i in range (len(listeCiblesAccessibles)):
        j=listeDistances[i][1]
        listeCiblesTriees.append(listeCiblesAccessibles[j])
    return listeCiblesTriees

def choixDeLaCible (xb,yb,listeDesCibles,theta): #Cette donction retourne les deux prochaines cibles du bateau
    listeDesCiblesCP = cp.copy(listeDesCibles)
    listeCiblesAccessibles=listeCiblesAcessibles(listeDesCiblesCP,xb,yb,theta)
    listeCiblesTriees = triCiblesDistance (xb,yb,listeCiblesAccessibles)


    for i in range (len(listeCiblesTriees)):
        listeDesCibles2 = []
        xpotentiel=listeCiblesTriees[i][0]
        ypotentiel =listeCiblesTriees [i][1]
        listeDesCibles2=listeDesCiblesCP
        listeDesCibles2.remove(listeCiblesTriees[i])
        listeCiblesAccessibles2=listeCiblesAcessibles(listeDesCibles2,xpotentiel,ypotentiel,theta)
        listeCiblesTriees2 = triCiblesDistance (xpotentiel,ypotentiel,listeCiblesAccessibles2)
        if listeCiblesTriees2!=[]:
            #print(1)
            return listeCiblesTriees[i]
        #print(2)
        return listeCiblesTriees[0]



## Simulation de l'asservissement et de la navigation
# COORDONNEES A ENTRER ICI__________________________________________________________________
xIni = 0
yIni = 0
thetaIni = 0
listeDesCibles = [[0,900],[-300,500],[0,-700],[500,0],[-300,0],[-800,-400]]
listeDesCibles1 = [[0,200]]
#____________________________________________________________________________________________


first = 0 #compteur pour vérifier le premier passage dans animate, ne pas modifier
c=10 #compteur pour la visibiliité des tournants (mettre à -1 pour enlever)
v = -1 #compteur pour le coup de vent (mettre à -1 pour enlever)


def coordBateauIntegration(xBateau , yBateau , thetaBateau, xCible, yCible):
    global c
    thetaCommande, puissanceHeliceAbsolue = asservissement(xCible , yCible , xBateau , yBateau , thetaBateau)
    #print('thetaCommande : '+str(thetaCommande))
    hasard = (100 + 10*((np.random.rand(1)[0] - 1/2)*2) )/100 #10% de hasard
    if abs(thetaCommande)>grandAngle:
        thetaReel = (thetaBateau + thetaCommande/2)%360
        vitesseReelle = puissanceHeliceAbsolue*(0.02)*hasard
    else :
        vitesseReelle = puissanceHeliceAbsolue*(0.02)*hasard
        if thetaCommande>0:
            thetaReel = (thetaBateau + grandAngle + (thetaCommande-grandAngle)/5)%360
        else :
            thetaReel = (thetaBateau - grandAngle + (thetaCommande+grandAngle)/5)%360

    #thetaReel = (thetaBateau + thetaCommande/2)%360
    #vitesseReelle = puissanceHeliceAbsolue*(0.01)*hasard #en m/s
    if c >0:
        vitesseReelle = 10
        c = c-1
    x = xBateau + cos(thetaBateau*(pi/180))*vitesseReelle*deltaT
    y = yBateau + sin(thetaBateau*(pi/180))*vitesseReelle*deltaT
    message = 'Le bateau est en (' + str(int(x)) + ',' + str(int(y)) + '). Angle : ' + str(int(thetaReel%360)) + '°. Vitesse : ' + str(round(vitesseReelle,2)) + 'm/s.'
    print(message)
    return x, y, thetaReel

##Trace animee de la simulation
fig, ax = plt.subplots()
x, y = [xIni],[yIni]
theta = [thetaIni]
plt.scatter(xIni, yIni, s=100, c="green")
for cible in listeDesCibles:
    plt.scatter(cible[0], cible[1], s=100, c="red")
sc = ax.scatter(x,y,s=5, c="black")
plt.xlim(-1000,1000)
plt.ylim(-1000,1000)


listeActuelle = []
def animate(i):
    global first
    global v
    global c
    global xCible
    global yCible
    global listeActuelle

    if first == 0: #initialisation
        first = first-1
        #listeActuelle = cp.copy(listeDesCibles)
        for e in listeDesCibles :
            listeActuelle.append(e)
        # nextFirst = choixDeLaCible(xIni,yIni,listeActuelle,thetaIni)
        # xCible = nextFirst[0]
        # yCible = nextFirst[1]

    xOld, yOld, thetaOld = x[-1], y[-1], theta[-1] #passage au point de trajectoire suivant
    next = choixDeLaCible(xOld,yOld,listeActuelle,thetaOld)
    if next == None :
        # if listeActuelle == []:
        #     ani.event_source.stop()
        xCible = listeActuelle[0][0]
        yCible = listeActuelle[0][1]
    else :
        xCible = next[0]
        yCible = next[1]
    l =coordBateauIntegration(xOld, yOld, thetaOld,xCible, yCible)
    x.append(l[0])
    y.append(l[1])
    theta.append(l[2])

    v = v-1
    if v == 0: #coup de vent
        x.append(x[-1]-100)
        y.append(y[-1]+150)
        theta.append((theta[-1]+90)%360)
        c = 10
        print("Coup de vent !")

    sc.set_offsets(np.c_[x,y])

    if sqrt((xCible-l[0])**2 + (yCible-l[1])**2)<rayonArrivee: #passage au point suivant
        if listeActuelle.count([xCible,yCible]) != 0:
            listeActuelle.remove([xCible,yCible])
        # if listeActuelle == []:
        #     ani.event_source.stop()
        c = 10
    if listeActuelle == []:
        ani.event_source.stop()

# def stopAni():
#     global listeActuelle
#     global first
#     if first == 0:
#         yield 0
#     else :
#         i=1
#         while listeActuelle != []:
#             i+=1
#             yield i

ani = matplotlib.animation.FuncAnimation(fig, animate, frames=1, interval=(deltaT*1), repeat=True)
plt.show()
























