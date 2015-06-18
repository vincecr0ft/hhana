# Imports
from ROOT import *
import numpy as np
import math
import myoptobs_2ndorder_noflav as HAWK
import myweightandoptobs as HAWK_new


try:
    import lhapdf
except ImportError:
    'print No lha pdf'
else:
    lhapdf.initPDFSet('CT10.LHgrid')

#----------------------
def sortJets(jet1,jet2,jet3):
    if jet1.Pt()<10e-7:
        jet1=TLorentzVector(0,0,0,0)
    if jet2.Pt()<10e-7:
        jet2=TLorentzVector(0,0,0,0)
    if jet3.Pt()<10e-7:
        jet3=TLorentzVector(0,0,0,0)
    if jet1.Pt()>jet2.Pt():
        lead=jet1
        sublead=jet2
    else:
        lead=jet2
        sublead=jet1
    if jet3.Pt()>lead.Pt():
        third=sublead
        sublead=lead
        lead=jet3
    elif jet3.Pt()>sublead.Pt():
        third=sublead
        sublead=jet3
    else:
        third=jet3
    return (lead,sublead,third)

def sortPartons(jet1,jet2,jet3,f1,f2,f3):
    if jet1.Pt()<10e-7:
        jet1=TLorentzVector(0,0,0,0)
    if jet2.Pt()<10e-7:
        jet2=TLorentzVector(0,0,0,0)
    if jet3.Pt()<10e-7:
        jet3=TLorentzVector(0,0,0,0)
    if f3 > 9 or f3 == None or f3 == 0.:
        if f3==21: 
            f3=0
        if jet1.Pt()>jet2.Pt():
            lead=[jet1,f1]
            sublead=[jet2,f2]
        else:
            lead=[jet2,f2]
            sublead=[jet1,f1]
        third=[jet3,f3]
    elif f2 > 9 or f2 == None or f2 == 0.:
        if f2==21: 
            f2=0
        if jet1.Pt()>jet3.Pt():
            lead=[jet1,f1]
            sublead=[jet3,f3]
        else:
            lead=[jet3,f3]
            sublead=[jet1,f1]
        third=[jet2,f2]
    elif f1 > 9 or f1 == None or f1 == 0.:
        if f1==21: 
            f1=0
        if jet2.Pt()>jet3.Pt():
            lead=[jet2,f2]
            sublead=[jet3,f3]
        else:
            lead=[jet3,f3]
            sublead=[jet2,f2]
        third=[jet1,f1]
    else:
        if jet1.Pt()>jet2.Pt():
            lead=[jet1,f1]
            sublead=[jet2,f2]
        else:
            lead=[jet2,f2]
            sublead=[jet1,f1]
        if jet3.Pt()>lead[0].Pt():
            third=sublead
            sublead=lead
            lead=[jet3,f3]
        elif jet3.Pt()>sublead[0].Pt():
            third=sublead
            sublead=[jet3,f3]
        else:
            third=[jet3,f3]
    return (lead,sublead,third)

def getX1(ecm,jet1,jet2,jet3,higgs):
    vec=jet1+jet2+jet3+higgs
    return (vec.M() * exp( vec.Rapidity() ) / ecm)

def getX2(ecm,jet1,jet2,jet3,higgs):
    vec=jet1+jet2+jet3+higgs
    return (vec.M() * exp( -vec.Rapidity() ) / ecm)

def getLHAPDF(Q,x):
    lhapdf.initPDF(0)
    pdfvals = np.array([])
    for f in range(-6,7):
        pdfvals = np.append(pdfvals,lhapdf.xfx(x,Q,f))
    return pdfvals

def getOptimalObservable(jet1,jet2,jet3,higgs):

    lead,sublead,third=sortJets(jet1,jet2,jet3)
    if lead == None or sublead == None or higgs == None: 
        print 'sort error'
        return None

    pjet1 = np.array([lead.E(),lead.Px(),lead.Py(),lead.Pz()])
    pjet2 = np.array([sublead.E(),sublead.Px(),sublead.Py(),sublead.Pz()])

    phiggs = np.array([higgs.E(),higgs.Px(),higgs.Py(),higgs.Pz()])

    ecm = 8000.
    x1 = getX1(ecm,jet1,jet2,jet3,higgs)
    x2 = getX2(ecm,jet1,jet2,jet3,higgs)


    if math.isnan(x1) or math.isnan(x2):
        print 'bjorken x error'
        return None
    
    if x1 == None or x2 == None:
        print 'bjorken x error'
        return None
    
    optobs = HAWK.optobs(ecm,x1,x2,pjet1,pjet2,phiggs)
    if math.isnan(optobs[0]):
        print 'optobs for: ecm %f, x1 %f, x2 %f, pjet1[0] %f, pjet2[0] %f, higgs %f'%(ecm,x1,x2,pjet1[0],pjet2[0],phiggs[0])
    return (optobs[0])

def getNewOptimalObservable(jet1,jet2,jet3,higgs):

    lead,sublead,third=sortJets(jet1,jet2,jet3)
    if lead == None or sublead == None or higgs == None: 
        print 'sort error'
        return None

    pjet1 = np.array([lead.E(),lead.Px(),lead.Py(),lead.Pz()])
    pjet2 = np.array([sublead.E(),sublead.Px(),sublead.Py(),sublead.Pz()])
    phiggs = np.array([higgs.E(),higgs.Px(),higgs.Py(),higgs.Pz()])


    if third == None or third.Pt()<10e-7:
        outCount = 2
    else:
        outCount = 3

    ecm = 8000.
    x1 = getX1(ecm,lead,sublead,third,higgs)
    x2 = getX2(ecm,lead,sublead,third,higgs)


    if math.isnan(x1) or math.isnan(x2):
        print 'bjorken x error'
        return None
    
    if x1 == None or x2 == None:
        print 'bjorken x error'
        return None

    pdfvals1 = getLHAPDF(125,x1)
    pdfvals2 = getLHAPDF(125,x2)

    return HAWK_new.optobs(ecm,higgs.M(),x1,x2,pdfvals1,pdfvals2,pjet1,pjet2,phiggs)[0]

def getDeltaPhi(jet1,jet2,jet3):
    lead, sublead, third = sortJets(jet1,jet2,jet3)
    if lead.Pt()<10e-7 or sublead.Pt()<10e-7 :return None

    if lead.Pz() > 0. and sublead.Pz() > 0.:
        dphi = None
    elif lead.Pz() < 0. and sublead. Pz() < 0.:
        dphi = None
    elif lead.Pz() > 0.:
        dphi = lead.DeltaPhi(sublead)
    else:
        dphi = sublead.DeltaPhi(lead)
    return dphi

def getOptimalWeight(d_tilde,jet1,jet2,jet3,higgs):

    lead,sublead,third=sortJets(jet1,jet2,jet3)
    if lead == None or sublead == None or higgs == None: return None

    pjet1 = np.array([lead.E(),lead.Px(),lead.Py(),lead.Pz()])
    pjet2 = np.array([sublead.E(),sublead.Px(),sublead.Py(),sublead.Pz()])

    phiggs = np.array([higgs.E(),higgs.Px(),higgs.Py(),higgs.Pz()])

    ecm = 8000.

    x1 = getX1(ecm,lead,sublead,third,higgs)
    x2 = getX2(ecm,lead,sublead,third,higgs)

    if x1 == None or x2 == None:
        return None
    weight = HAWK.reweight(ecm,1,1,0,0,d_tilde,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x1,x2,pjet1,pjet2,phiggs)
    if not math.isnan(weight):
        return weight
    else:
        return 0.

def getNewOptimalWeight(d_tilde,jet1,jet2,jet3,higgs,
                        partin1_f,partin2_f,
                        parton1_f,parton2_f,parton3_f):

    lead,sublead,third=sortPartons(jet1,jet2,jet3,parton1_f,parton2_f,parton3_f)

    if lead[0].Pt() < 10e-7 or sublead[0].Pt() < 10e-7 or higgs.M() <10e-7: 
        print 'no lead sublead or higgs'
        print lead[0].Pt(),' ',sublead[0].Pt(),' ',higgs.M()
        return None

    pjet1 = np.array([lead[0].E(),lead[0].Px(),lead[0].Py(),lead[0].Pz()])
    pjet2 = np.array([sublead[0].E(),sublead[0].Px(),sublead[0].Py(),sublead[0].Pz()])
    phiggs = np.array([higgs.E(),higgs.Px(),higgs.Py(),higgs.Pz()])

    if third[0].Pt()<10e-7:
        pjet3 = np.array([0,0,0,0])
        outCount = 2
        third=[TLorentzVector(0,0,0,0),0]
    else:
        pjet3 = np.array([third[0].E(),third[0].Px(),third[0].Py(),third[0].Pz()])
        outCount = 3
    ecm = 8000.
    x1 = getX1(ecm,jet1,jet2,jet3,higgs)
    x2 = getX2(ecm,jet1,jet2,jet3,higgs)
    if x1 == None or x2 == None:
        print 'bjorken x fail'
        return None
    if third[1] == None:
        third[1] = 0.

    if partin1_f==21:
        partin1_f=0
    if partin2_f==21:
        partin2_f=0
    ModelSwitch = 1 # 0 for pure CP-odd, 1 for CP-mixed
    rsmin = 1
    din = 0
    dbin = 0
    lambdahvvin = -1
    try :
        weight = HAWK_new.reweight(ecm,higgs.M(),ModelSwitch,rsmin,din,dbin,d_tilde,
                           d_tilde,lambdahvvin,0,0,0,0,0,0,0,0,0,0,0,0,outCount,
                           partin1_f,partin2_f,
                           lead[1],sublead[1],third[1],
                           x1,x2,pjet1,pjet2,pjet3,phiggs)[0]
    except:
        print 'weight failed'
        print 'ecm %f,ModelSwitch %f,rsmin %f,din%f ,dbin %f,d_tilde %f'%(ecm,ModelSwitch,rsmin,din,dbin,d_tilde)
        print 'lamda %f, outCount %f'%(lambdahvvin,outCount)
        print 'partin1 f %f, partin2 %f'%(partin1_f,partin2_f)
        print 'lead[1] %f,sublead[1] %f,third[1] %f,'%(lead[1],sublead[1],third[1])
        print 'x1 %f, x2 %f, pjet1 %f, pjet2 %f, pjet3 %f, phiggs %f'%(x1,x2,pjet1[0],pjet2[0],pjet3[0],phiggs[0])
#    print 'hawkweight of ',weight
    if weight<0.:
        print 'negative weight'
        print 'ecm %f,ModelSwitch %f,rsmin %f,din%f ,dbin %f,d_tilde %f'%(ecm,ModelSwitch,rsmin,din,dbin,d_tilde)
        print 'lamda %f, outCount %f'%(lambdahvvin,outCount)
        print 'partin1 f %f, partin2 %f'%(partin1_f,partin2_f)
        print 'lead[1] %f,sublead[1] %f,third[1] %f,'%(lead[1],sublead[1],third[1])
        print 'x1 %f, x2 %f, pjet1 %f, pjet2 %f, pjet3 %f, phiggs %f'%(x1,x2,pjet1[0],pjet2[0],pjet3[0],phiggs[0])
    if not math.isnan(weight):
#        print 'incoming 1 f %f, incoming 2 %f'%(partin1_f,partin2_f)
#        print 'lead flav %f,sublead flav %f,third flav %f,'%(lead[1],sublead[1],third[1])


        return weight
    else:
        print 'nanananananananana nan weight'
        return None


def buildEvent(
        jet1_pt,
        jet2_pt,
        jet3_pt,
        jet1_eta,
        jet2_eta,
        jet3_eta,
        jet1_phi,
        jet2_phi,
        jet3_phi,
        tau1_pt,
        tau1_eta,
        tau1_phi,
        tau1_m,
        numJets):

    jet1=TLorentzVector(0.,0.,0.,0.)
    jet2=TLorentzVector(0.,0.,0.,0.)
    jet3=TLorentzVector(0.,0.,0.,0.)
    higgs=TLorentzVector(0.,0.,0.,0.)
    jet1.SetPtEtaPhiM(jet1_pt*0.001,jet1_eta,jet1_phi,0.)
    jet2.SetPtEtaPhiM(jet2_pt*0.001,jet2_eta,jet2_phi,0.)
    jet3.SetPtEtaPhiM(jet3_pt*0.001,jet3_eta,jet3_phi,0.)
    higgs.SetPtEtaPhiM(tau1_pt,tau1_eta,tau1_phi,tau1_m)
    
    nJets=0
    if jet1.Pt()>10e-7:nJets+=1
    if jet2.Pt()>10e-7:nJets+=1
    if jet3.Pt()>10e-7:nJets+=1

    if numJets<2 or nJets<2:
        return None
    elif higgs.M()<10e-7:
        print 'mmc m too low ',higgs.M()
        return None
    else:
        return getOptimalObservable(jet1,jet2,jet3,higgs)

def buildNewEvent(
        jet1_pt,
        jet2_pt,
        jet3_pt,
        jet1_eta,
        jet2_eta,
        jet3_eta,
        jet1_phi,
        jet2_phi,
        jet3_phi,
        tau1_pt,
        tau1_eta,
        tau1_phi,
        tau1_m,
        numJets):

    jet1=TLorentzVector(0.,0.,0.,0.)
    jet2=TLorentzVector(0.,0.,0.,0.)
    jet3=TLorentzVector(0.,0.,0.,0.)
    higgs=TLorentzVector(0.,0.,0.,0.)
    jet1.SetPtEtaPhiM(jet1_pt*0.001,jet1_eta,jet1_phi,0.)
    jet2.SetPtEtaPhiM(jet2_pt*0.001,jet2_eta,jet2_phi,0.)
    jet3.SetPtEtaPhiM(jet3_pt*0.001,jet3_eta,jet3_phi,0.)
    higgs.SetPtEtaPhiM(tau1_pt,tau1_eta,tau1_phi,tau1_m)
    
    nJets=0
    if jet1.Pt()>10e-7:nJets+=1
    if jet2.Pt()>10e-7:nJets+=1
    if jet3.Pt()>10e-7:nJets+=1

    if numJets<2 or nJets<2:
        return None
    elif higgs.M()<10e-7:
        print 'mmc m too low ',higgs.M()
        return None
    else:
        return getNewOptimalObservable(jet1,jet2,jet3,higgs)


def buildPhiEvent(
        jet1_pt,
        jet2_pt,
        jet3_pt,
        jet1_eta,
        jet2_eta,
        jet3_eta,
        jet1_phi,
        jet2_phi,
        jet3_phi,
        numJets):

    jet1=TLorentzVector(0.,0.,0.,0.)
    jet2=TLorentzVector(0.,0.,0.,0.)
    jet3=TLorentzVector(0.,0.,0.,0.)
    jet1.SetPtEtaPhiM(jet1_pt*0.001,jet1_eta,jet1_phi,0.)
    jet2.SetPtEtaPhiM(jet2_pt*0.001,jet2_eta,jet2_phi,0.)
    jet3.SetPtEtaPhiM(jet3_pt*0.001,jet3_eta,jet3_phi,0.)
    nJets=0
    if jet1.Pt()>10e-7:nJets+=1
    if jet2.Pt()>10e-7:nJets+=1
    if jet3.Pt()>10e-7:nJets+=1

    if numJets<2 or nJets<2:
        return None
    else:
        dphi = getDeltaPhi(jet1,jet2,jet3)
        return dphi

def buildEventWeight(d_tilde,
        jet1_pt,
        jet2_pt,
        jet3_pt,
        jet1_eta,
        jet2_eta,
        jet3_eta,
        jet1_phi,
        jet2_phi,
        jet3_phi,
        jet1_m,
        jet2_m,
        jet3_m,
        tau1_pt,
        tau1_eta,
        tau1_phi,
        tau1_m,
        numJets,
        weight):

    jet1=TLorentzVector(0.,0.,0.,0.)
    jet2=TLorentzVector(0.,0.,0.,0.)
    jet3=TLorentzVector(0.,0.,0.,0.)
    higgs=TLorentzVector(0.,0.,0.,0.)
    jet1.SetPtEtaPhiM(jet1_pt*0.001,jet1_eta,jet1_phi,0.)
    jet2.SetPtEtaPhiM(jet2_pt*0.001,jet2_eta,jet2_phi,0.)
    jet3.SetPtEtaPhiM(jet3_pt*0.001,jet3_eta,jet3_phi,0.)
    higgs.SetPtEtaPhiM(tau1_pt,tau1_eta,tau1_phi,tau1_m)


    if numJets<2 or higgs.M()<10e-7: 
        return weight
    else:
        try: hawkweight = getOptimalWeight(d_tilde,jet1,jet2,jet3,higgs)
        except: 
            print "failed trying weight "
            return weight
    if math.isnan(hawkweight):
        print "failed weight ",hawkweight
        return weight
    else:
        return weight*hawkweight

def buildNewEventWeight(d_tilde,
        jet1_pt,
        jet2_pt,
        jet3_pt,
        jet1_eta,
        jet2_eta,
        jet3_eta,
        jet1_phi,
        jet2_phi,
        jet3_phi,
        tau1_pt,
        tau1_eta,
        tau1_phi,
        tau1_m,
        partin1_f,
        partin2_f,
        parton1_f,
        parton2_f,
        parton3_f,
        numJets,
        weight):
    if tau1_m < 10e-7:
        return weight

    jet1=TLorentzVector(0.,0.,0.,0.)
    jet2=TLorentzVector(0.,0.,0.,0.)
    jet3=TLorentzVector(0.,0.,0.,0.)
    higgs=TLorentzVector(0.,0.,0.,0.)
    jet1.SetPtEtaPhiM(jet1_pt*0.001,jet1_eta,jet1_phi,0.)
    jet2.SetPtEtaPhiM(jet2_pt*0.001,jet2_eta,jet2_phi,0.)
    jet3.SetPtEtaPhiM(jet3_pt*0.001,jet3_eta,jet3_phi,0.)
    higgs.SetPtEtaPhiM(tau1_pt*0.001,tau1_eta,tau1_phi,tau1_m*0.001)

    nJets=0
    if jet1.Pt()>10e-7:nJets+=1
    if jet2.Pt()>10e-7:nJets+=1
    if jet3.Pt()>10e-7:nJets+=1

    if numJets<2 or nJets<2:
        return weight
    try:
        hawkweight = getNewOptimalWeight(d_tilde,jet1,jet2,jet3,higgs,
                                           partin1_f,partin2_f,
                                           parton1_f,parton2_f,parton3_f)
    except Exception as e: 
        print "failed trying weight "
        print 'exception: ',e
        print '(d_tilde %f,jet1 %f,jet2 %f,jet3 %f,higgs %f)'%(d_tilde,jet1.Pt(),jet2.Pt(),jet3.Pt(),higgs.M())
        print '(partin1_f %i,partin2_f %i,parton1_f %i,parton2_f %i,parton3_f %i)'%(partin1_f,partin2_f, parton1_f,parton2_f,parton3_f)
        return weight
#
    if hawkweight == None:
        print 'no hawkweight'
        return weight
    elif math.isnan(hawkweight):
        print "nan weight ",hawkweight
        return weight
    elif hawkweight < 0.:
        print "negative weight!"
        return weight
    else:
        return weight*hawkweight


observize=np.vectorize(buildEvent)
new_observize=np.vectorize(buildNewEvent)
dphize=np.vectorize(buildPhiEvent)
weight_vize=np.vectorize(buildEventWeight)
new_weight_vize=np.vectorize(buildNewEventWeight)
