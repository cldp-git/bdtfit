#include <TString.h>
class ConfigParams {
  public :
	int m_nBinsData;
	int m_markerStyle;
	int m_font;
    int m_axisfont;
	int m_canvWidth;
	int m_canvHeight;
	double m_legendMargin;
	double m_legendTextSize;
	double m_legendPos[4];
    double m_titlesize;
	
	//Begin definition of shared physics parameters
	double m_cbnp_bd;
	double m_cbnp_bs;
	double m_cbnm_bd;
	double m_cbnm_bs;
	double m_cbap_bd;
	double m_cbap_bs;
	double m_cbam_bd;
	double m_cbam_bs;
	double m_signal_sigmap_bs;
	double m_signal_sigmap_bd;
	double m_signal_sigmam_bs;
	double m_signal_sigmam_bd;
	double m_fracm_bs;
	double m_fracm_bd;
	
	double m_mcwidthfactor_signal;

	double m_massshift_2012_MC;
    double m_lbmassshift_2013_MC;
    double m_bdmassshift_2013_MC;
	double m_mcwidthfactor_2012_MC;
	
	double m_bsbd_mdiff;
	double m_bubd_mdiff;
    double m_bdlb_mdiff; 
    bool m_fixcombfromws;
	bool m_floathorns;

	bool m_usegaus;
	double m_gausfrac;
        double m_siggauswidth;

	bool m_dotoys;
	bool m_toyvariation_extragausgenerated;
	float m_toyvariation_extragausgenerated_fraction,m_toyvariation_extragausgenerated_width;
	
	bool m_cutonfdchi2;
	double m_fdchi2cut;

	double m_fracm_all;	
	double m_cbnp_all;
	double m_cbnm_all;
	double m_cbap_all;
	double m_cbam_all;

	double fixedChlessReflForDsDs;
	double fixedLcReflForDsDs;
	double fixedDReflForDsDs;
    double fixedChlessReflForLcDs;
        bool fixChlessReflForDsDs;
	bool fixLcReflForDsDs;
	bool fixDReflForDsDs;

	bool fixChlessForOthers;

	TString m_suffix;

	ConfigParams(){
		m_nBinsData = 110;
		m_markerStyle = 24;
		m_font = 42;
        m_axisfont = 42;
		m_canvWidth = 750;
		m_canvHeight = 500;
		m_legendMargin = 0.25;
		m_legendTextSize = .06;//0.045;
		m_legendPos[0] = 0.65;
		m_legendPos[1] = 0.475;
		m_legendPos[2] = 0.90;
		m_legendPos[3] = 0.90;
		m_titlesize = .06;
		//Begin definition of shared physics parameters
		m_cbnp_bd = 10.;
		m_cbnp_bs = 10.;
		m_cbnm_bd = 10.;
		m_cbnm_bs = 10.;
		m_cbap_bd = -1.64;//-1.51;
		m_cbap_bs = -.97;//-1.42;
		m_cbam_bd = 1.23;//0.69;
		m_cbam_bs = 0.68;//1.04;
		m_signal_sigmap_bd = 7.18;
		m_signal_sigmap_bs = 6.60;
		m_signal_sigmam_bd = 7.40;
		m_signal_sigmam_bs = 8.06;
		m_fracm_bs	   = 0.31;
		m_fracm_bd	   = 0.13;	

		m_mcwidthfactor_signal = 1.2;

		m_massshift_2012_MC = 1.9;
        m_lbmassshift_2013_MC = .8;
        m_bdmassshift_2013_MC = .7;
		m_mcwidthfactor_2012_MC = 1.1;
		
		m_bsbd_mdiff = 87.3;
		m_bubd_mdiff = 0.33;
		m_bdlb_mdiff = 339.2;//339.2;
		m_fixcombfromws = true;
		m_floathorns = true;
		
		m_fracm_all = m_fracm_bs;
		m_cbnm_all  = m_cbnm_bs;
		m_cbnp_all  = m_cbnp_bs;
		m_cbam_all  = m_cbam_bs;
		m_cbap_all  = m_cbap_bs;

		m_suffix = "_DTF";

		m_cutonfdchi2 = true;
		m_fdchi2cut = 2.;

		m_dotoys = false;
		m_toyvariation_extragausgenerated = true;
		m_toyvariation_extragausgenerated_fraction = .16;
		m_toyvariation_extragausgenerated_width = 14.4;

		m_usegaus = false;
		m_gausfrac = .16;
 	        m_siggauswidth = 14.4;

		fixedChlessReflForDsDs = 0.01;
        fixedChlessReflForLcDs = 0.034;
		fixedDReflForDsDs = 35.;
		fixedLcReflForDsDs = 15.;

		fixDReflForDsDs = true;
		fixChlessReflForDsDs = true;
		fixLcReflForDsDs = true;

		fixChlessForOthers = true;

	};
	inline virtual ~ConfigParams() {};

  
  private:

  ClassDef(ConfigParams,1)
};

