#include <string>

#include "TFile.h"
#include "TNtuple.h"
#include "TNtupleD.h"
#include "TSystem.h"
#include "TPad.h"
#include "TStyle.h"

#include "lcfiplus.h"
#include "process.h"
#include "DNNProvider2.h"
#include "VertexSelector.h"
#include "algoEtc.h"
#include "VertexFinderSuehara.h"
#include "VertexFitterSimple.h"

#include <utility>
#include <algorithm>

using namespace lcfiplus;

namespace lcfiplus {

void DNNProvider2::init(Parameters* param) {
  Algorithm::init(param);
  string filename = param->get("DNNProvider2.FileName",string("DNNProvider2.root"));
  _jetname = param->get("DNNProvider2.JetCollectionName",string("RefinedJets"));
  string privtx = param->get("PrimaryVertexCollectionName",string("PrimaryVertex"));
  Event::Instance()->setDefaultPrimaryVertex(privtx.c_str());

  _mcIsB = param->get("DNNProvider2.MC_IsB",int(0));
  _mcIsC = param->get("DNNProvider2.MC_IsC",int(0));
  _mcIsU = param->get("DNNProvider2.MC_IsU",int(0));
  _mcIsD = param->get("DNNProvider2.MC_IsD",int(0));
  _mcIsS = param->get("DNNProvider2.MC_IsS",int(0));
  _mcIsG = param->get("DNNProvider2.MC_IsG",int(0)); 

  _jets = 0;

  if(_mcIsB + _mcIsC + _mcIsU + _mcIsD + _mcIsS + _mcIsG != 1){
    cout << "Label parameter is wrong! Put 1 on one of MC_Is(BCQ)." << endl;
  }

  cout << filename << endl;

  _file = new TFile(filename.c_str(),"RECREATE");
  _ntp = new TTree("tree","tree");

  DNNData & d = _data;
  // jet variables
  _ntp->Branch("jet_px",&d.jet_px,"jet_px/F");
  _ntp->Branch("jet_py",&d.jet_py,"jet_py/F");
  _ntp->Branch("jet_pz",&d.jet_pz,"jet_pz/F");
  _ntp->Branch("jet_mass",&d.jet_mass,"jet_mass/F");
  _ntp->Branch("jet_ntracks",&d.jet_ntracks,"jet_ntracks/I");
  _ntp->Branch("jet_nneutrals",&d.jet_nneutrals,"jet_nneutrals/I");

  _ntp->Branch("jet_phi",&d.jet_phi,"jet_phi/F");
  _ntp->Branch("jet_theta",&d.jet_theta,"jet_theta/F");

  // for charges
  // particle kinematics
  _ntp->Branch("pfcand_px",&d.px);
  _ntp->Branch("pfcand_py",&d.py);
  _ntp->Branch("pfcand_pz",&d.pz);
  _ntp->Branch("pfcand_e",&d.e);
  _ntp->Branch("pfcand_efrac",&d.efrac);
  _ntp->Branch("pfcand_erel_log",&d.erel_log);
  _ntp->Branch("pfcand_thetarel",&d.dtheta);
  _ntp->Branch("pfcand_phirel",&d.dphi);
  _ntp->Branch("pfcand_thetarel_ilc",&d.dtheta_ilc);
  _ntp->Branch("pfcand_phirel_ilc",&d.dphi_ilc);

  // track errors
  _ntp->Branch("pfcand_dptdpt",&d.cov_omega);
  _ntp->Branch("pfcand_detadeta",&d.cov_tanlambda);
  _ntp->Branch("pfcand_dphidphi",&d.cov_phi);
  _ntp->Branch("pfcand_dxydxy",&d.cov_d0);
  _ntp->Branch("pfcand_dzdz",&d.cov_z0);
  _ntp->Branch("pfcand_dxydz",&d.cov_d0_z0);
  _ntp->Branch("pfcand_dphidxy",&d.cov_d0_phi);
  _ntp->Branch("pfcand_dlambdadz",&d.cov_z0_tanlambda);
  _ntp->Branch("pfcand_dxyc",&d.cov_d0_omega);
  _ntp->Branch("pfcand_dxyctgtheta",&d.cov_d0_tanlambda);
  _ntp->Branch("pfcand_phic",&d.cov_phi_omega);
  _ntp->Branch("pfcand_phidz",&d.cov_z0_phi);
  _ntp->Branch("pfcand_phictgtheta",&d.cov_phi_tanlambda);
  _ntp->Branch("pfcand_cdz",&d.cov_z0_omega);
  _ntp->Branch("pfcand_cctgtheta",&d.cov_omega_tanlambda);

  // particle displacements
  _ntp->Branch("d0",&d.d0);
  _ntp->Branch("d0sig",&d.d0sig);
  _ntp->Branch("z0",&d.z0);
  _ntp->Branch("z0sig",&d.z0sig);
  _ntp->Branch("ip3d",&d.ip3d);
  _ntp->Branch("ip3dsig",&d.ip3dsig);

  _ntp->Branch("dEdx",&d.dEdx);

  _ntp->Branch("pfcand_dxy",&d.dxy);
  _ntp->Branch("pfcand_dz",&d.dz);
  _ntp->Branch("pfcand_btagSip2dVal",&d.ip2d_fcc);
  _ntp->Branch("pfcand_btagSip2dSig",&d.ip2dsig_fcc);
  _ntp->Branch("pfcand_btagSip3dVal",&d.ip3d_fcc);
  _ntp->Branch("pfcand_btagSip3dSig",&d.ip3dsig_fcc);
  _ntp->Branch("pfcand_btagJetDistVal",&d.jetdist_fcc);
  _ntp->Branch("pfcand_btagJetDistSig",&d.jetdistsig_fcc);

  // particle ID
  _ntp->Branch("pfcand_charge",&d.charge);
  _ntp->Branch("pfcand_isMu",&d.ismuon);
  _ntp->Branch("pfcand_isEl",&d.iselectron);
  _ntp->Branch("pfcand_isGamma",&d.isphoton);
  _ntp->Branch("pfcand_isChargedHad",&d.ischargedhadron);
  _ntp->Branch("pfcand_isNeutralHad",&d.isneutralhadron);
  _ntp->Branch("pfcand_type",&d.pdg_pfa);
  _ntp->Branch("pfcand_mcpid",&d.mcpid);
  _ntp->Branch("pfcand_mcp_pdg",&d.mcp_pdg);
  _ntp->Branch("pfcand_Ktype",&d.K_pdg_pfa); //20240201
  _ntp->Branch("pfcand_isPion",&d.ispion);
  _ntp->Branch("pfcand_isKaon",&d.iskaon);
  _ntp->Branch("pfcand_isProton",&d.isproton);
  // _ntp->Branch("pfcand_isKaon0", &d.iskaon0);

  // for PI3 (20240203)
  _ntp->Branch("pfcand_proton_K",&d.proton_K);
  _ntp->Branch("pfcand_pion_K",&d.pion_K);
  _ntp->Branch("pfcand_proton_Klike",&d.proton_Klike);
  _ntp->Branch("pfcand_pion_Klike",&d.pion_Klike);
  _ntp->Branch("pfcand_dEdxEl",&d.electron_dEdxdistance);
  _ntp->Branch("pfcand_dEdxMu",&d.muon_dEdxdistance);
  _ntp->Branch("pfcand_dEdxPion",&d.pion_dEdxdistance);
  _ntp->Branch("pfcand_dEdxKaon",&d.kaon_dEdxdistance);
  _ntp->Branch("pfcand_dEdxProton",&d.proton_dEdxdistance);
  _ntp->Branch("pfcand_iselectronlike",&d.iselectronlike);
  _ntp->Branch("pfcand_ismuonlike",&d.ismuonlike);
  _ntp->Branch("pfcand_ispionlike",&d.ispionlike);
  _ntp->Branch("pfcand_iskaonlike",&d.iskaonlike);
  _ntp->Branch("pfcand_isprotonlike",&d.isprotonlike);




  // for neutrals
  // particle kinematics
  _ntp->Branch("neu_pfcand_px",&d.neu_px);
  _ntp->Branch("neu_pfcand_py",&d.neu_py);
  _ntp->Branch("neu_pfcand_pz",&d.neu_pz);
  _ntp->Branch("neu_pfcand_e",&d.neu_e);
  _ntp->Branch("neu_pfcand_efrac",&d.neu_efrac);
  _ntp->Branch("neu_pfcand_erel_log",&d.neu_erel_log);
  _ntp->Branch("neu_pfcand_thetarel",&d.neu_dtheta);
  _ntp->Branch("neu_pfcand_phirel",&d.neu_dphi);
  _ntp->Branch("neu_pfcand_thetarel_ilc",&d.neu_dtheta_ilc);
  _ntp->Branch("neu_pfcand_phirel_ilc",&d.neu_dphi_ilc);

  // track errors
  _ntp->Branch("neu_pfcand_dptdpt",&d.neu_cov_omega);
  _ntp->Branch("neu_pfcand_detadeta",&d.neu_cov_tanlambda);
  _ntp->Branch("neu_pfcand_dphidphi",&d.neu_cov_phi);
  _ntp->Branch("neu_pfcand_dxydxy",&d.neu_cov_d0);
  _ntp->Branch("neu_pfcand_dzdz",&d.neu_cov_z0);
  _ntp->Branch("neu_pfcand_dxydz",&d.neu_cov_d0_z0);
  _ntp->Branch("neu_pfcand_dphidxy",&d.neu_cov_d0_phi);
  _ntp->Branch("neu_pfcand_dlambdadz",&d.neu_cov_z0_tanlambda);
  _ntp->Branch("neu_pfcand_dxyc",&d.neu_cov_d0_omega);
  _ntp->Branch("neu_pfcand_dxyctgtheta",&d.neu_cov_d0_tanlambda);
  _ntp->Branch("neu_pfcand_phic",&d.neu_cov_phi_omega);
  _ntp->Branch("neu_pfcand_phidz",&d.neu_cov_z0_phi);
  _ntp->Branch("neu_pfcand_phictgtheta",&d.neu_cov_phi_tanlambda);
  _ntp->Branch("neu_pfcand_cdz",&d.neu_cov_z0_omega);
  _ntp->Branch("neu_pfcand_cctgtheta",&d.neu_cov_omega_tanlambda);

  // particle displacements
  _ntp->Branch("neu_d0",&d.neu_d0);
  _ntp->Branch("neu_d0sig",&d.neu_d0sig);
  _ntp->Branch("neu_z0",&d.neu_z0);
  _ntp->Branch("neu_z0sig",&d.neu_z0sig);
  _ntp->Branch("neu_ip3d",&d.neu_ip3d);
  _ntp->Branch("neu_ip3dsig",&d.neu_ip3dsig);

  _ntp->Branch("neu_pfcand_dxy",&d.neu_dxy);
  _ntp->Branch("neu_pfcand_dz",&d.neu_dz);
  _ntp->Branch("neu_pfcand_btagSip2dVal",&d.neu_ip2d_fcc);
  _ntp->Branch("neu_pfcand_btagSip2dSig",&d.neu_ip2dsig_fcc);
  _ntp->Branch("neu_pfcand_btagSip3dVal",&d.neu_ip3d_fcc);
  _ntp->Branch("neu_pfcand_btagSip3dSig",&d.neu_ip3dsig_fcc);
  _ntp->Branch("neu_pfcand_btagJetDistVal",&d.neu_jetdist_fcc);
  _ntp->Branch("neu_pfcand_btagJetDistSig",&d.neu_jetdistsig_fcc);

  // particle ID
  _ntp->Branch("neu_pfcand_charge",&d.neu_charge);
  _ntp->Branch("neu_pfcand_isMu",&d.neu_ismuon);
  _ntp->Branch("neu_pfcand_isEl",&d.neu_iselectron);
  _ntp->Branch("neu_pfcand_isGamma",&d.neu_isphoton);
  _ntp->Branch("neu_pfcand_isChargedHad",&d.neu_ischargedhadron);
  _ntp->Branch("neu_pfcand_isNeutralHad",&d.neu_isneutralhadron);
  _ntp->Branch("neu_pfcand_type",&d.neu_pdg_pfa);
  _ntp->Branch("neu_pfcand_mcpid",&d.neu_mcpid);
  _ntp->Branch("pfcand_mcp_pdg",&d.mcp_pdg);
  _ntp->Branch("neu_pfcand_Ktype",&d.neu_K_pdg_pfa); //20240201
  _ntp->Branch("neu_pfcand_isPion",&d.neu_ispion);
  _ntp->Branch("neu_pfcand_isKaon",&d.neu_iskaon);
  _ntp->Branch("neu_pfcand_isProton",&d.neu_isproton);
  // _ntp->Branch("neu_pfcand_isKaon0", &d.neu_iskaon0);


  // label
  _ntp->Branch("mc_b",&d.mc_b,"mc_b/I");
  _ntp->Branch("mc_c",&d.mc_c,"mc_c/I");
  _ntp->Branch("mc_u",&d.mc_u,"mc_u/I");
  _ntp->Branch("mc_d",&d.mc_d,"mc_d/I");
  _ntp->Branch("mc_s",&d.mc_s,"mc_s/I");
  _ntp->Branch("mc_g",&d.mc_g,"mc_g/I");
  _ntp->Branch("mc_q",&d.mc_q,"mc_q/I"); // usdg

}

void DNNProvider2::process() {
  if (!_jets) {
    Event::Instance()->Get(_jetname.c_str(), _jets);
  }

  const JetVec& jets = *_jets;
  const Vertex *privtx = Event::Instance()->getPrimaryVertex();

  DNNData &d = _data;

  for (unsigned int j=0; j < jets.size(); ++j) {
    const Jet* jet = jets[j];

    memset(&_data,0,sizeof(_data));

    d.jet_px = jet->Px();
    d.jet_py = jet->Py();
    d.jet_pz = jet->Pz();
    d.jet_e = jet->E();
    d.jet_mass = jet->M();
    TrackVec &tracks = jet->getAllTracks();
    NeutralVec &neutrals = jet->getNeutrals();
    // MCParticleVec &mcps = jet->getMCParticle();

    d.jet_theta = jet->Theta();
    d.jet_phi = jet->Phi();

    d.jet_ntracks = tracks.size();
    d.jet_nneutrals = neutrals.size();

    //float jet_theta = jet->Theta();
    //float jet_phi = jet->Phi();

    d.mc_b = _mcIsB;
    d.mc_c = _mcIsC;
    d.mc_u = _mcIsU;
    d.mc_d = _mcIsD;
    d.mc_s = _mcIsS;
    d.mc_g = _mcIsG;
    d.mc_q = _mcIsU || _mcIsD || _mcIsS || _mcIsG;

    // probably order of tracks/netural does not matter...
    int nall = d.jet_ntracks + d.jet_nneutrals;
    int ntra = d.jet_ntracks;
    int nneu = d.jet_nneutrals;
    if (ntra==0) continue;
    if (nneu==0) continue;
    //for charges
    d.px.resize(ntra);
    d.py.resize(ntra);
    d.pz.resize(ntra);
    d.e.resize(ntra);
    d.efrac.resize(ntra);
    d.erel_log.resize(ntra);
    d.dtheta.resize(ntra);
    d.dphi.resize(ntra);
    d.dtheta_ilc.resize(ntra);
    d.dphi_ilc.resize(ntra);

    d.cov_d0.resize(ntra);
    d.cov_z0.resize(ntra);
    d.cov_phi.resize(ntra);
    d.cov_omega.resize(ntra);
    d.cov_tanlambda.resize(ntra);

    d.cov_d0_z0.resize(ntra);
    d.cov_d0_phi.resize(ntra);
    d.cov_d0_omega.resize(ntra);
    d.cov_d0_tanlambda.resize(ntra);

    d.cov_z0_phi.resize(ntra);
    d.cov_z0_omega.resize(ntra);
    d.cov_z0_tanlambda.resize(ntra);

    d.cov_phi_omega.resize(ntra);
    d.cov_phi_tanlambda.resize(ntra);
    d.cov_omega_tanlambda.resize(ntra);

    d.d0.resize(ntra);
    d.d0sig.resize(ntra);
    d.z0.resize(ntra);
    d.z0sig.resize(ntra);
    d.ip3d.resize(ntra);
    d.ip3dsig.resize(ntra);

    d.dEdx.resize(ntra);

    d.dxy.resize(ntra);
    d.dz.resize(ntra);
    d.ip2d_fcc.resize(ntra);
    d.ip2dsig_fcc.resize(ntra);
    d.ip3d_fcc.resize(ntra);
    d.ip3dsig_fcc.resize(ntra);
    d.jetdist_fcc.resize(ntra);
    d.jetdistsig_fcc.resize(ntra);

    d.charge.resize(ntra);
    d.ismuon.resize(ntra);
    d.iselectron.resize(ntra);
    d.isphoton.resize(ntra);
    d.ischargedhadron.resize(ntra);
    d.isneutralhadron.resize(ntra);
    d.ispion.resize(ntra);
    d.iskaon.resize(ntra);
    // d.iskaon0.resize(ntra);
    d.isproton.resize(ntra);
    d.pdg_pfa.resize(ntra);
    d.mcpid.resize(ntra);
    d.mcp_pdg.resize(ntra);
    d.K_pdg_pfa.resize(ntra);

    d.proton_K.resize(ntra);
    d.pion_K.resize(ntra);
    d.proton_Klike.resize(ntra);
    d.pion_Klike.resize(ntra);

    d.electron_dEdxdistance.resize(ntra);
    d.muon_dEdxdistance.resize(ntra);
    d.pion_dEdxdistance.resize(ntra);
    d.kaon_dEdxdistance.resize(ntra);
    d.proton_dEdxdistance.resize(ntra);

    d.iselectronlike.resize(ntra);
    d.ismuonlike.resize(ntra);
    d.ispionlike.resize(ntra);
    d.iskaonlike.resize(ntra);
    d.isprotonlike.resize(ntra);



    // for neutrals
    d.neu_px.resize(nneu);
    d.neu_py.resize(nneu);
    d.neu_pz.resize(nneu);
    d.neu_e.resize(nneu);
    d.neu_efrac.resize(nneu);
    d.neu_erel_log.resize(nneu);
    d.neu_dtheta.resize(nneu);
    d.neu_dphi.resize(nneu);
    d.neu_dtheta_ilc.resize(nneu);
    d.neu_dphi_ilc.resize(nneu);

    d.neu_cov_d0.resize(nneu);
    d.neu_cov_z0.resize(nneu);
    d.neu_cov_phi.resize(nneu);
    d.neu_cov_omega.resize(nneu);
    d.neu_cov_tanlambda.resize(nneu);

    d.neu_cov_d0_z0.resize(nneu);
    d.neu_cov_d0_phi.resize(nneu);
    d.neu_cov_d0_omega.resize(nneu);
    d.neu_cov_d0_tanlambda.resize(nneu);

    d.neu_cov_z0_phi.resize(nneu);
    d.neu_cov_z0_omega.resize(nneu);
    d.neu_cov_z0_tanlambda.resize(nneu);

    d.neu_cov_phi_omega.resize(nneu);
    d.neu_cov_phi_tanlambda.resize(nneu);
    d.neu_cov_omega_tanlambda.resize(nneu);

    d.neu_d0.resize(nneu);
    d.neu_d0sig.resize(nneu);
    d.neu_z0.resize(nneu);
    d.neu_z0sig.resize(nneu);
    d.neu_ip3d.resize(nneu);
    d.neu_ip3dsig.resize(nneu);

    d.neu_dxy.resize(nneu);
    d.neu_dz.resize(nneu);
    d.neu_ip2d_fcc.resize(nneu);
    d.neu_ip2dsig_fcc.resize(nneu);
    d.neu_ip3d_fcc.resize(nneu);
    d.neu_ip3dsig_fcc.resize(nneu);
    d.neu_jetdist_fcc.resize(nneu);
    d.neu_jetdistsig_fcc.resize(nneu);

    d.neu_charge.resize(nneu);
    d.neu_ismuon.resize(nneu);
    d.neu_iselectron.resize(nneu);
    d.neu_isphoton.resize(nneu);
    d.neu_ischargedhadron.resize(nneu);
    d.neu_isneutralhadron.resize(nneu);
    d.neu_ispion.resize(nneu);
    d.neu_iskaon.resize(nneu);
    // d.neu_iskaon0.resize(nneu);
    d.neu_isproton.resize(nneu);
    d.neu_pdg_pfa.resize(nneu);
    d.neu_mcpid.resize(nneu);
    d.neu_mcp_pdg.resize(nneu);
    d.neu_K_pdg_pfa.resize(nneu);


    vector<std::pair<float, int> > order_tr;
    order_tr.resize(ntra);
    vector<std::pair<float, int> > order_n;
    order_n.resize(nneu);

    int i;

    for(i=0;i<d.jet_ntracks;i++){
      if(ntra==0){
        cout << "passedT" << endl;
        continue;
      }
      const Track *tr = tracks[i];
      order_tr[i] = std::pair<float, int>(tr->E(), i);
    }
    for(i=0;i<nneu;i++){
      if(nneu==0){
        cout << "passedN" << endl;
        continue;
      }
      const Neutral *neu = neutrals[i];
      order_n[i] = std::pair<float, int>(neu->E(), i);
    }

    // sort by energy
    //   std::sort(order.begin(), order.end(), [](std::pair<float, int>a, std::pair<float, int> b){
	  // return a.first > b.first;
    //     });



    std::sort(order_tr.begin(), order_tr.end(), [](std::pair<float, int>a, std::pair<float, int> b){
	    return a.first > b.first;
      });
    std::sort(order_n.begin(), order_n.end(), [](std::pair<float, int>a, std::pair<float, int> b){
	    return a.first > b.first;
      });



    for(i=0;i<ntra;i++){
      // if(order[i].second >= d.jet_ntracks) continue;
      //cout << i << " " << order[i].second << " " << d.jet_ntracks << " " << nall << endl;
      if(ntra==0){
        cout << "passedT" << endl;
        continue;
      }

      const Track *tr = tracks[order_tr[i].second];
      // const MCParticle *mcparticles = mcps[order_tr[i].second];
      // const Track *tr = tracks[i];
      d.px[i] = tr->Px();
      d.py[i] = tr->Py();
      d.pz[i] = tr->Pz();
      d.e[i] = tr->E();
      d.efrac[i] = tr->E() / jet->E();
      d.erel_log[i] = log10(d.efrac[i]);
      d.dtheta_ilc[i] = tr->Theta() - jet->Theta();
      d.dphi_ilc[i] = tr->Phi() - jet->Phi();
      if(d.dphi_ilc[i] < -TMath::Pi())d.dphi_ilc[i] += TMath::Pi() * 2;
      if(d.dphi_ilc[i] > TMath::Pi())d.dphi_ilc[i] -= TMath::Pi() * 2;
      calc_thetaphi(jet->Vect(), tr->Vect(), d.dtheta[i], d.dphi[i]);

      // track covmatrix
      d.cov_d0[i] = tr->getCovMatrix()[tpar::d0d0];
      d.cov_z0[i] = tr->getCovMatrix()[tpar::z0z0];
      d.cov_phi[i] = tr->getCovMatrix()[tpar::phph];
      d.cov_omega[i] = tr->getCovMatrix()[tpar::omom];
      d.cov_tanlambda[i] = tr->getCovMatrix()[tpar::tdtd];

      d.cov_d0_z0[i] = tr->getCovMatrix()[tpar::d0z0];
      d.cov_d0_phi[i] = tr->getCovMatrix()[tpar::d0ph];
      d.cov_d0_omega[i] = tr->getCovMatrix()[tpar::d0om];
      d.cov_d0_tanlambda[i] = tr->getCovMatrix()[tpar::d0td];

      d.cov_z0_phi[i] = tr->getCovMatrix()[tpar::z0ph];
      d.cov_z0_omega[i] = tr->getCovMatrix()[tpar::z0om];
      d.cov_z0_tanlambda[i] = tr->getCovMatrix()[tpar::z0td];

      d.cov_phi_omega[i] = tr->getCovMatrix()[tpar::phom];
      d.cov_phi_tanlambda[i] = tr->getCovMatrix()[tpar::phtd];
      d.cov_omega_tanlambda[i] = tr->getCovMatrix()[tpar::omtd];

      d.d0[i] = tr->getD0();
      d.d0sig[i] = tr->getD0() / sqrt(tr->getCovMatrix()[tpar::cov::d0d0]);
      d.z0[i] = tr->getZ0();
      d.z0sig[i] = tr->getZ0() / sqrt(tr->getCovMatrix()[tpar::cov::z0z0]);

      d.ip3d[i] = sqrt(tr->getD0() * tr->getD0() + tr->getZ0() * tr->getZ0());
      d.ip3dsig[i] = d.ip3d[i] / sqrt(tr->getCovMatrix()[tpar::cov::d0d0] + tr->getCovMatrix()[tpar::cov::z0z0] + 2 * tr->getCovMatrix()[tpar::cov::d0z0]);

      d.dEdx[i] = tr->getdEdx();
      
      d.dEdx[i] = tr->getdEdx();
      
      d.dxy[i] = calc_dxy(tr->getD0(), tr->getZ0(), tr->getPhi(), tr->Vect(), privtx->getPos(), tr->getCharge());
      d.dz[i] = calc_dz(tr->getD0(), tr->getZ0(), tr->getPhi(), tr->Vect(), privtx->getPos(), tr->getCharge());
      d.ip2d_fcc[i] = calc_sip2d(tr->getD0(), tr->getPhi(), jet->Px(), jet->Py());
      d.ip2dsig_fcc[i] = d.ip2d_fcc[i] / sqrt(d.cov_d0[i]);
      d.ip3d_fcc[i] = calc_sip3d(tr->getD0(), tr->getZ0(), tr->getPhi(), jet->Vect());
      d.ip3dsig_fcc[i] = d.ip3d_fcc[i] / sqrt(d.cov_d0[i] + d.cov_z0[i]);
      d.jetdist_fcc[i] = calc_jetDist(tr->getD0(), tr->getZ0(), tr->getPhi(), tr->Vect(), jet->Vect());
      d.jetdistsig_fcc[i] = d.jetdist_fcc[i] / sqrt(d.cov_d0[i] + d.cov_z0[i]);

      d.charge[i] = tr->getCharge();

      // get ParticleID 0 or 1
      // d.ischargedhadron[i] = !(d.ismuon[i] || d.iselectron[i]);
      // d.isneutralhadron[i] = 0.0;
      // std::vector<float> charged(5);
      // charged.clear();
      // d.ismuon[i] = 0;
      // d.iselectron[i] = 0;
      // d.iskaon[i] = 0;
      // d.ispion[i] = 0;
      // d.isproton[i] = 0;
      // d.isphoton[i] = 0;
      // charged.push_back(tr->getParticleIDProbability("muonProbability"));
      // charged.push_back (tr->getParticleIDProbability("electronProbability"));
      // charged.push_back(tr->getParticleIDProbability("kaonProbability"));
      // charged.push_back(tr->getParticleIDProbability("pionProbability"));
      // charged.push_back(tr->getParticleIDProbability("protonProbability"));

      // std::vector<float>::iterator iter = std::max_element(charged.begin(), charged.end());
      // int index = std::distance(charged.begin(), iter);
      // if (index==0) d.ismuon[i] = 1;
      // if (index==1) d.iselectron[i] = 1;
      // if (index==2) d.iskaon[i] = 1;
      // if (index==3) d.ispion[i] = 1;
      // if (index==4) d.isproton[i] = 1;
      // d.pdg_pfa[i] = tr->getPDG();


      // tracing LCFIPlus default
      d.isphoton[i] = 0;
      d.ismuon[i] = tr->getParticleIDProbability("muonProbability");

      d.iselectron[i] = tr->getParticleIDProbability("electronProbability");
      d.isphoton[i] = 0.0;
      d.ischargedhadron[i] = !(d.ismuon[i] || d.iselectron[i]);
      d.isneutralhadron[i] = 0.0;
      d.ispion[i] = tr->getParticleIDProbability("pionProbability");
      d.iskaon[i] = tr->getParticleIDProbability("kaonProbability");
      d.isproton[i] = tr->getParticleIDProbability("protonProbability");
      d.pdg_pfa[i] = tr->getPDG();

      if(tr->getMcp()){
	d.mcpid[i] = (tr->getMcp())->getId();
	d.mcp_pdg[i] = (tr->getMcp())->getPDG();
      }

      d.K_pdg_pfa[i] = 0;
      if (d.pdg_pfa[i]==321) d.K_pdg_pfa[i] = 1;
      if (d.pdg_pfa[i]==310) d.K_pdg_pfa[i] = 1;
      d.pion_K[i] = d.ispion[i] - d.iskaon[i];
      d.proton_K[i] = d.isproton[i] - d.iskaon[i];
      d.pion_Klike[i] = tr->getParticleIDProbability("pionLikelihood") - tr->getParticleIDProbability("kaonLikelihood");
      d.proton_Klike[i] = tr->getParticleIDProbability("protonLikelihood") - tr->getParticleIDProbability("kaonLikelihood");
      
      d.electron_dEdxdistance[i] = tr->getParticleIDProbability("electron_dEdxdistance");
      d.muon_dEdxdistance[i] = tr->getParticleIDProbability("muon_dEdxdist ance");
      d.pion_dEdxdistance[i] = tr->getParticleIDProbability("pion_dEdxdistance");
      d.kaon_dEdxdistance[i] = tr->getParticleIDProbability("kaon_dEdxdistance");
      d.proton_dEdxdistance[i] = tr->getParticleIDProbability("proton_dEdxdistance");

      d.iselectronlike[i] = tr->getParticleIDProbability("electronLikelihood");
      d.ismuonlike[i] = tr->getParticleIDProbability("muonLikelihood");
      d.ispionlike[i] = tr->getParticleIDProbability("pionLikelihood");
      d.iskaonlike[i] = tr->getParticleIDProbability("kaonLikelihood");
      d.isprotonlike[i] = tr->getParticleIDProbability("protonLikelihood");

      // PDG ID

      // d.pdg_pfa[i] = tr->getPDG();
      // int pidnum = tr->getPDG();  
      // d.ismuon[i] = 0;
      // d.iselectron[i] = 0;
      // d.iskaon[i] = 0;
      // d.ispion[i] = 0;
      // d.isproton[i] = 0;
      // d.isphoton[i] = 0;
      // d.isneutralhadron[i] = 0;
      // d.iskaon0[i] = 0;
      // d.ischargedhadron[i] = !(d.ismuon[i] || d.iselectron[i]);
      // if (pidnum==11){
      //   d.iselectron[i] = 1;
      //   d.ischargedhadron[i] = !(d.ismuon[i] || d.iselectron[i]);
      // }
      // if (pidnum==13){
      //   d.ismuon[i] = 1;
      //   d.ischargedhadron[i] = !(d.ismuon[i] || d.iselectron[i]);
      // }
      // if (pidnum==2112){
      //   d.isneutralhadron[i] = 1;
      // }
      // if (pidnum==22){
      //   d.isphoton[i] = 1;
      // }
      // if (pidnum==211){
      //   d.ispion[i] = 1;
      // }
      // if (pidnum==321){
      //   d.iskaon[i] = 1;
      // }
      // if (pidnum==310){
      //   d.isneutralhadron[i] = 1;
      //   d.iskaon0[i] = 1;
      // }
      // if (pidnum==2212){
      //   d.isproton[i] = 1;
      // }


     // original version
      // d.ismuon[i] = algoEtc::SimpleSecMuonFinder(tr, 5,5,5, -0.1, 0.2, 0.8, 1.5, 4, 0.5, privtx);
      // d.iselectron[i] = algoEtc::SimpleSecElectronFinder(tr, 5,5,5,5,0.98,0.9, 1.15, privtx);
      // d.isphoton[i] = 0;
      // d.ischargedhadron[i] = !(d.ismuon[i] || d.iselectron[i]);
      // d.isneutralhadron[i] = 0;
      // d.pdg_pfa[i] = tr->getPDG();

    }
    // for neutrals
    for(i=0;i<nneu;i++){
      // if(order[i].second < d.jet_ntracks) continue;
      if(nneu==0){
        cout << "passedN" << endl;
        continue;
      }
      // const Neutral *neu = neutrals[i];
      const Neutral *neu = neutrals[order_n[i].second];
      // const MCParticle *mcparticles = mcps[order_n[i].second];
      d.neu_px[i] = neu->Px();
      d.neu_py[i] = neu->Py();
      d.neu_pz[i] = neu->Pz();
      d.neu_e[i] = neu->E();
      d.neu_efrac[i] = neu->E() / jet->E();
      d.neu_erel_log[i] = log10(d.neu_efrac[i]);
      d.neu_dtheta_ilc[i] = neu->Theta() - jet->Theta();
      d.neu_dphi_ilc[i] = neu->Phi() - jet->Phi();
      if(d.neu_dphi_ilc[i] < -TMath::Pi())d.neu_dphi_ilc[i] += TMath::Pi() * 2;
      if(d.neu_dphi_ilc[i] > TMath::Pi())d.neu_dphi_ilc[i] -= TMath::Pi() * 2;
      calc_thetaphi(jet->Vect(), neu->Vect(), d.neu_dtheta[i], d.neu_dphi[i]);

      d.neu_cov_d0[i] = -9;
      d.neu_cov_z0[i] = -9;
      d.neu_cov_phi[i] = -9;
      d.neu_cov_omega[i] = -9;
      d.neu_cov_tanlambda[i] = -9;

      d.neu_cov_d0_z0[i] = -9;
      d.neu_cov_d0_phi[i] = -9;
      d.neu_cov_d0_omega[i] = -9;
      d.neu_cov_d0_tanlambda[i] = -9;

      d.neu_cov_z0_phi[i] = -9;
      d.neu_cov_z0_omega[i] = -9;
      d.neu_cov_z0_tanlambda[i] = -9;

      d.neu_cov_phi_omega[i] = -9;
      d.neu_cov_phi_tanlambda[i] = -9;
      d.neu_cov_omega_tanlambda[i] = -9;

      d.neu_dxy[i] = -9;
      d.neu_dz[i] = -9;
      d.neu_ip2d_fcc[i] = -9;
      d.neu_ip2dsig_fcc[i] = -9;
      d.neu_ip3d_fcc[i] = -9;
      d.neu_ip3dsig_fcc[i] = -9;
      d.neu_jetdist_fcc[i] = -9;
      d.neu_jetdistsig_fcc[i] = -9;

      d.neu_charge[i] = 0;
      d.neu_ismuon[i] = 0;
      d.neu_iselectron[i] = 0;
      // simple photon finder
      double ecaldep = neu->getCaloEdep()[tpar::ecal];
      double hcaldep = neu->getCaloEdep()[tpar::hcal];
      d.neu_isphoton[i] = (ecaldep / (ecaldep + hcaldep) > 0.98);
      d.neu_ischargedhadron[i] = 0;
      d.neu_isneutralhadron[i] = !d.isphoton[i];
      d.neu_ispion[i] = 0;
      d.neu_iskaon[i] = 0;
      d.neu_isproton[i] = 0;
      // d.neu_iskaon0[i] = 0;
      d.neu_pdg_pfa[i] = neu->getPDG();
      if(neu->getMcp()){
	d.neu_mcpid[i] = (neu->getMcp())->getId();
	d.neu_mcp_pdg[i] = (neu->getMcp())->getPDG();
      }
      d.neu_K_pdg_pfa[i] = 0;
      if (d.neu_pdg_pfa[i]==321) d.neu_K_pdg_pfa[i] = 1;
      if (d.neu_pdg_pfa[i]==310) d.neu_K_pdg_pfa[i] = 1;
    
      // d.neu_ismuon[i] = neu->getParticleIDProbability("muonProbability");
      // d.neu_iselectron[i] = neu->getParticleIDProbability("electronProbability");
      // d.neu_ispion[i] = neu->getParticleIDProbability("pionProbability");
      // d.neu_iskaon[i] = neu->getParticleIDProbability("kaonProbability");
      // d.neu_isproton[i] = neu->getParticleIDProbability("protonProbability");


      // PDG ID
      // d.pdg_pfa[i] = neu->getPDG();
      // int pidnum = neu->getPDG();  
      // d.ismuon[i] = 0;
      // d.iselectron[i] = 0;
      // d.iskaon[i] = 0;
      // d.ispion[i] = 0;
      // d.isproton[i] = 0;
      // d.isphoton[i] = 0;
      // d.isneutralhadron[i] = 0;
      // d.iskaon0[i] = 0;
      // d.ischargedhadron[i] = !(d.ismuon[i] || d.iselectron[i]);
      // if (pidnum==11){
      //   d.iselectron[i] = 1;
      //   d.ischargedhadron[i] = !(d.ismuon[i] || d.iselectron[i]);
      // }
      // if (pidnum==13){
      //   d.ismuon[i] = 1;
      //   d.ischargedhadron[i] = !(d.ismuon[i] || d.iselectron[i]);
      // }
      // if (pidnum==2112){
      //   d.isneutralhadron[i] = 1;
      // }
      // if (pidnum==22){
      //   d.isphoton[i] = 1;
      // }
      // if (pidnum==211){
      //   d.ispion[i] = 1;
      // }
      // if (pidnum==321){
      //   d.iskaon[i] = 1;
      // }
      // if (pidnum==310){
      //   d.isneutralhadron[i] = 1;
      //   d.iskaon0[i] = 1;
      // }
      // if (pidnum==2212){
      //   d.isproton[i] = 1;
      // }
    

    }

    _ntp->Fill();
  }
}

void DNNProvider2::end() {
  _file->Write();
  _file->Close();
}

// copied from FCCANalyses/analyzers/dataframe/src/ReconstructedParticle2Track.cc

float DNNProvider2::calc_dxy(float D0_wrt0, float Z0_wrt0, float phi0_wrt0, TVector3 p, TVector3 privtx, int charge){
  double Bz = 3.5;
  double cSpeed = 2.99792458e8 * 1e-9;

  TVector3 X( - D0_wrt0 * TMath::Sin(phi0_wrt0) , D0_wrt0 * TMath::Cos(phi0_wrt0) , Z0_wrt0);
  TVector3 x = X - privtx;
  //std::cout<<"vertex: "<<V.Vect().X()<<", "<<V.Vect().Y()<<", "<<V.Vect().Z()<<", "<<std::endl;

  double a = - charge * Bz * cSpeed;
  double pt = p.Pt();
  double r2 = x(0) * x(0) + x(1) * x(1);
  double cross = x(0) * p(1) - x(1) * p(0);
  double D=-9;
  if (pt * pt - 2 * a * cross + a * a * r2 > 0) {
    double T = TMath::Sqrt(pt * pt - 2 * a * cross + a * a * r2);
    if (pt < 10.0) D = (T - pt) / a;
    else D = (-2 * cross + a * r2) / (T + pt);
  }
  return D;
}

float DNNProvider2::calc_dz(float D0_wrt0, float Z0_wrt0, float phi0_wrt0, TVector3 p, TVector3 privtx, int charge){
  double Bz = 3.5;
  double cSpeed = 2.99792458e8 * 1e-9;

  TVector3 X( - D0_wrt0 * TMath::Sin(phi0_wrt0) , D0_wrt0 * TMath::Cos(phi0_wrt0) , Z0_wrt0);
  TVector3 x = X - privtx;

  double a = - charge * Bz * cSpeed;
  double pt = p.Pt();
  double C = a/(2 * pt);
  double r2 = x(0) * x(0) + x(1) * x(1);
  double cross = x(0) * p(1) - x(1) * p(0);
  double T = TMath::Sqrt(pt * pt - 2 * a * cross + a * a * r2);
  double D;
  if (pt < 10.0) D = (T - pt) / a;
  else D = (-2 * cross + a * r2) / (T + pt);
  double B = C * TMath::Sqrt(TMath::Max(r2 - D * D, 0.0) / (1 + 2 * C * D));
  if ( TMath::Abs(B) > 1.) B = TMath::Sign(1, B);
  double st = TMath::ASin(B) / C;
  double ct = p(2) / pt;
  double z0;
  double dot = x(0) * p(0) + x(1) * p(1);
  if (dot > 0.0) z0 = x(2) - ct * st;
  else z0 = x(2) + ct * st;

  return z0;
}

float DNNProvider2::calc_sip2d(float D0, float phi0, float jetpx, float jetpy)
{
  TVector2 p(jetpx, jetpy);

  TVector2 d0(-D0 * TMath::Sin(phi0), D0 * TMath::Cos(phi0));
  return TMath::Sign(1, d0 * p) * fabs(D0);
}

float DNNProvider2::calc_sip3d(float D0, float Z0, float phi0, TVector3 p_jet)
{
  TVector3 d(-D0 * TMath::Sin(phi0), D0 * TMath::Cos(phi0), Z0);
  return TMath::Sign(1, d * p_jet) * fabs(sqrt(D0 * D0 + Z0 * Z0));
}

float DNNProvider2::calc_jetDist(float D0, float Z0, float phi0, TVector3 p_ct, TVector3 p_jet)
{
  TVector3 d(-D0 * TMath::Sin(phi0), D0 * TMath::Cos(phi0), Z0);
  TVector3 r_jet(0.0, 0.0, 0.0);
  TVector3 n = p_ct.Cross(p_jet).Unit(); // What if they are parallel?
  return n.Dot(d - r_jet);
}

void DNNProvider2::calc_thetaphi(TVector3 jet, TVector3 part, float &theta, float &phi){
  part.RotateZ(-jet.Phi());
  part.RotateY(-jet.Theta());

  theta = part.Theta();
  phi = part.Phi();  
}

}

