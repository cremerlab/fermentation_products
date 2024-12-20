#basic caclulations
import os
import json
import numpy as np


def estimate_carbsforbacteria(total_carbon,fiber,sugar,fiberdigestion=0.75,starchpassage=0.15):
    return fiberdigestion*fiber+starchpassage*(total_carbon-fiber-sugar)
#function to transform fecal wet weight into dry weight
#this is based on data points from many different studies (as derived below)
#there is a correlation between dry and weg weight but the functional form is hard to estimate
#here diferent fits are used ('calculationmodes')

#    hadzadiet["energyBE05"],fermc,hadzadiet['fermentationprodBE05'],fermg,fermg2,cc,order=energycalc(hadzadiet["carbLI05"],scenario='reference',calctype='from_carbs')
#rewrite having only one output

def energycalc_array(inputv,scenario='reference',calctype='none',dict_yielddata="example.json"):
    outputv_energy=inputv.copy()
    outputv_ferm=inputv.copy()
    for il in range(0,inputv.shape[0]):
        outputv_energy[il],c1,outputv_ferm[il],c3,c4,c5,c6,c7=energycalc(inputv[il],scenario=scenario,calctype=calctype,dict_yielddata=dict_yielddata)
    return outputv_energy,outputv_ferm




def cal_energy_FP(excretion,order):
    energyc=0
    enthalpy_sub=['glucose','maltose','acetate','butyrate','formate','lactate','propionate','succinate']
    enthalpy=[0.68,1.36,0.21,0.52,0.,.33,0.37,0.36] 
    iS = -1
    #print(order)
    #print(excretion)
    for sub in order:
        iS=iS+1
        energyc=energyc+excretion[iS]*enthalpy[enthalpy_sub.index(sub)]
    return energyc

#estimate amount of water produced, O2 consumed, CO2 released.

def cal_CO2_o2(total_carbon,fiber,protein,fat,fiberdigestion=0.75,starchpassage=0.15,proteinpassage=0.07):

    #input in grams
    
    #estimate carbohydrates
    MAC=fiberdigestion*fiber+starchpassage*(total_carbon-fiber-sugar)
    nonMAC=total_carbon-MAC
    carbohydrates_respiration=nonMAC #assume all CH not reaching large intestine get respirate
    #fats
    fat_respiration=fat #assume all fats are absorbed and respirated
    #proteins
    proteins_respiration=(1.-proteinpassage)*protein

    #calculate oxygen consumption (in g)
    o2_consumption=2.87*fat_respiration+1.42*proteins_respiration+1.0656*carbohydrates_respiration #double check, number for oxygen consumption for fat

    #calculate oxygen
    
    #calculate water production (in g)
    h2o_release=1.112*fat_respiration+0.45*proteins_respiration+0.60*carbohydrates_respiration
    
    #calculate CO2 production (in g)
    co2_release=2.75*fat_respiration+1.74*proteins_respiration+1.465*carbohydrates_respiration

    #urea production
    urea_release=0.3*proteins_respiration

    return [o2_consumption,h2o_release,co2_release,urea_release]
    #carbohydrates
    #Mass of CO₂ produced: 146.52 g
	#Mass of H₂O produced: 59.99 g

    #fats
    #In summary, the complete oxidation of 100 g of fat (using palmitic acid as a representative) results in approximately 275 g of CO₂ and 112 g of water. Since fats primarily consist of carbon, hydrogen, and oxygen, they do not produce nitrogenous waste, unlike proteins.

    #    So, the metabolism of 100 g of protein results in approximately 174 g of CO₂, 45 g of water, and 30 g of urea.
    
#function to transform fecal wet weight into dry weight
#this is based on data points from many different studies (as derived below)
#there is a correlation between dry and weg weight but the functional form is hard to estimate
#here diferent fits are used ('calculationmodes')


def energycalc(inputv,scenario='reference',calctype='from_carbs',dict_yielddata="example.json"): #amount of carbohydrates reaching gut, in g/day; fraction of B.theta in population (assumes rest is E. rectale)
    
    with open(os.path.join(dict_yielddata)) as f:
        yielddict = json.load(f)
    
    inverseyieldref=0.18015*yielddict["uptake"]
    
    #take input as fecal bacterial dry mass or caculate bacterial drymass
    if calctype=='from_carbs':
        bacdrymass=inputv/inverseyieldref
    elif calctype=='from_feces':
        bacdrymass=inputv
    else:
        error_calctype_not_found

    sublist=['acetate','butyrate','formate','lactate','propionate','succinate']
    order=sublist
    sublist_color=['#1b9e77','#66a61e','#a6761d','#e7298a','#d95f02','#7570b3']
    sublist_conversion_to_gram=np.array([60.052,88.11,46.03,90.08,74.079,118.09])
    sublist_catoms=np.array([2,4,1,3,3,4])
    #sublist energy
    #sublist_energy=[0.21,0.52,0.,.33,0.37,0.36]
    #enthalpielist=[0.21,0.37,0.36,0.33,0.52,0.061,0.327]#kcal/mmol
    
    excretion=[]
    for sub in sublist:
        excretion.append(yielddict[sub])
    excretion=np.array(excretion)
    
    #amount of fermentation products
    totfermlist=bacdrymass*excretion
    totfermlist_sum=np.sum(totfermlist)
    
    #convert to gram
    totfermlist_gram=np.multiply(totfermlist,sublist_conversion_to_gram)/1000.
    totfermlist_gram_sum=np.sum(totfermlist_gram)
    #total energy
    totenergy=cal_energy_FP(excretion,order)*bacdrymass
    #previously: yielddict["total_secretion_energy"]
    #print("double check enthalpies!!!")
    
    #total c atoms mmol
    totfem_carb=np.sum(np.multiply(totfermlist,sublist_catoms))
    #totfem_carb=np.sum(np.multiply(totfermlist,totfermlist_gram))
    
    
    #np.sum(np.multiply(totfermlist,sublist_catoms))
    #print(totenergy)
    #print(totfermlist)
    
    return totenergy,totfermlist.tolist(),totfermlist_sum,totfermlist_gram.tolist(),totfermlist_gram_sum,totfem_carb,order,bacdrymass #returns energy of fermentation products (in kcal) and amount of fermentation products (in mmol)


def estimate_carbsforbacteria(total_carbon,fiber,sugar,fiberdigestion=0.75,starchpassage=0.15):
    return fiberdigestion*fiber+starchpassage*(total_carbon-fiber-sugar)
 
#function to transform fecal wet weight into dry weight
#this is based on data points from many different studies (as derived below)
#there is a correlation between dry and weg weight but the functional form is hard to estimate
#here diferent fits are used ('calculationmodes')
def dryweight_fromwetweight(wetmassin,calculationmode=1,fraction=False):
    if calculationmode==0:
        #constant fraction
        curc=0.25
    elif calculationmode==1:
        #linear decrease
        parc=[4.63384822e-04, 6.85305059e-01]
        curc=1.-(wetmassin*parc[0]+parc[1])
    elif calculationmode==3:
        #sqrt behavior + linear
        parc=[0.49246012,  0.03076236, -0.00070195]
        curc= 1.-(parc[0] +parc[1]*np.sqrt(wetmassin)+parc[2]*wetmassin)
    else:
        error
    if fraction==True:
        return curc
    else:
        return curc*wetmassin