===========
XSteam
===========
XSteam provides (mostly) accurate steam and water properties from 0 - 1000 bar and from 0 - 2000 °C according to the standard IAPWS IF-97.
For accuracy of the functions in different regions see IF-97 Page 4
(www.iapws.org) http://www.iapws.org/relguide/IF97-Rev.pdf
Basic usage::
    #!/usr/bin/env python

    from XSteam import XSteam

    steamTable = XSteam(mksSystem = True)
    print steamTable.hL_p(220.0)


Requirements
=========

Tests require numpy, Demos require numpy and matplotlib

Install
=========
First, make shure that numpy is installed
....
To test is setup was successful, run pyXSteamDemo.py from the comandline.
There is still (v0.0) some Error in all Regions, so be warned that almost all Regions excede the Error Range.
Apart form accuracy errors, there should be no warnings.

Nomenclature
=========
All Functions follow the same naming schema:
        First the wanted property, then a underscore "_", then the wanted input properties
        Example: t_ph is temperature as a function of pressure and enthalpy. 
        For a list of valid functions se bellow:
        t    Temperature	(°C or °F)
        p    Pressure	(bar or psi)
        h    Enthalpy	(kJ/kg or btu/lb)
        v    Specific volume	(m3/kg or ft^3/lb)
        rho  Density	(kg/m3 or lb/ft^3)
        s    Specific entropy	(kJ/(kg °C) or btu/(lb °F))
        u    Specific internal energy	(kJ/kg or btu/lb)
        Cp   Specific isobaric heat capacity	(kJ/(kg °C) or btu/(lb °F))
        Cv   Specific isochoric heat capacity	(kJ/(kg °C) or btu/(lb °F))
        w    Speed of sound	(m/s or ft/s)
        my   Viscosity	(N s/m^2 or lbm/ft/hr)
        tc   Thermal Conductivity	(W/(m °C) or btu/(h ft °F))
        st   Surface Tension	(N/m or lb/ft)
        x    Vapour fraction
        vx   Vapour Volume Fraction

Usage
=========
    By using the mksSystem Parameter, you can tell XSteam witch Unit System to use.
		mksSystem = True  : m/kg/sec/°C/bar/W
		mksSystem = False : ft/lb/sec/°F/psi/btu

    to run all Test:
	>>> from XSteam.XSteamTester import XSteamTester
	>>> steamTableTester = XSteamTester(verbose = True)
	>>> steamTableTester.runTests()
		
	Available Functions:
		Temperature	
			tsat_p	Saturation temperature
			t_ph	Temperture as a function of pressure and enthalpy
			t_ps	Temperture as a function of pressure and entropy
			t_hs	Temperture as a function of enthalpy and entropy
		
		Pressure	
			psat_t	Saturation pressure
			p_hs	Pressure as a function of h and s. 
			p_hrho Pressure as a function of h and rho. Very unaccurate for solid water region since it's almost incompressible!
		
		Enthalpy	
			hV_p	Saturated vapour enthalpy
			hL_p	Saturated liquid enthalpy
			hV_t	Saturated vapour enthalpy
			hL_t	Saturated liquid enthalpy
			h_pt	Entalpy as a function of pressure and temperature.
			h_ps	Entalpy as a function of pressure and entropy.
			h_px	Entalpy as a function of pressure and vapour fraction
			h_prho	Entalpy as a function of pressure and density. Observe for low temperatures (liquid) this equation has 2 solutions. 
			h_tx	Entalpy as a function of temperature and vapour fraction
		
		Specific volume	
			vV_p	Saturated vapour volume
			vL_p	Saturated liquid volume
			vV_t	Saturated vapour volume
			vL_t	Saturated liquid volume
			v_pT	Specific volume as a function of pressure and temperature.
			v_ph	Specific volume as a function of pressure and enthalpy
			v_ps	Specific volume as a function of pressure and entropy.
		
		Density	
			rhoV_p	Saturated vapour density
			rhoL_p	Saturated liquid density
			rhoV_t	Saturated vapour density
			rhoL_t	Saturated liquid density
			rho_pt	Density as a function of pressure and temperature.
			rho_ph	Density as a function of pressure and enthalpy
			rho_ps	Density as a function of pressure and entropy.
		
		Specific entropy 	
			sV_p	Saturated vapour entropy
			sL_p	Saturated liquid entropy
			sV_t	Saturated vapour entropy
			sL_t	Saturated liquid entropy
			s_pt	Specific entropy as a function of pressure and temperature (Returns saturated vapour entalpy if mixture.)
			s_ph	Specific entropy as a function of pressure and enthalpy
		
		Specific internal energy 	
			uV_p	Saturated vapour internal energy
			uL_p	Saturated liquid internal energy
			uV_t	Saturated vapour internal energy
			uL_t	Saturated liquid internal energy
			u_pT	Specific internal energy as a function of pressure and temperature.
			u_ph	Specific internal energy as a function of pressure and enthalpy
			u_ps	Specific internal energy as a function of pressure and entropy.
		
		Specific isobaric heat capacity 	
			CpV_p	Saturated vapour heat capacity 
			CpL_p	Saturated liquid heat capacity 
			CpV_t	Saturated vapour heat capacity 
			CpL_t	Saturated liquid heat capacity 
			Cp_pt	Specific isobaric heat capacity as a function of pressure and temperature.
			Cp_ph	Specific isobaric heat capacity as a function of pressure and enthalpy
			Cp_ps	Specific isobaric heat capacity as a function of pressure and entropy.
		
		Specific isochoric heat capacity 	
			CvV_p	Saturated vapour isochoric heat capacity
			CvL_p	Saturated liquid isochoric heat capacity
			CvV_t	Saturated vapour isochoric heat capacity
			CvL_t	Saturated liquid isochoric heat capacity
			Cv_pt	Specific isochoric heat capacity as a function of pressure and temperature.
			Cv_ph	Specific isochoric heat capacity as a function of pressure and enthalpy
			Cv_ps	Specific isochoric heat capacity as a function of pressure and entropy
		
		Speed of sound 	
			wV_p	Saturated vapour speed of sound
			wL_p	Saturated liquid speed of sound
			wV_t	Saturated vapour speed of sound
			wL_t	Saturated liquid speed of sound
			w_pt	Speed of sound as a function of pressure and temperature
			w_ph	Speed of sound as a function of pressure and enthalpy
			w_ps	Speed of sound as a function of pressure and entropy
		
		Viscosity	
			Viscosity is not part of IAPWS Steam IF97. Equations from 
			"Revised Release on the IAPWS Formulation 1985 for the Viscosity of Ordinary Water Substance", 2003 are used.	
			Viscosity in the mixed region (4) is interpolated according to the density. This is not true since it will be two fases.	
			my_pt	Viscosity as a function of pressure and temperature.
			my_ph	Viscosity as a function of pressure and enthalpy
			my_ps	Viscosity as a function of pressure and entropy
		
		Thermal Conductivity	
			Revised release on the IAPS Formulation 1985 for the Thermal Conductivity of ordinary water substance (IAPWS 1998)	
			tcL_p	Saturated vapour thermal conductivity
			tcV_p	Saturated liquid thermal conductivity
			tcL_t	Saturated vapour thermal conductivity
			tcV_t	Saturated liquid thermal conductivity
			tc_pT	Thermal conductivity as a function of pressure and temperature
			tc_ph	Thermal conductivity as a function of pressure and enthalpy
			tc_hs	Thermal conductivity as a function of enthalpy and entropy
		
		Surface tension
			st_t	Surface tension for two phase water/steam as a function of T
			st_p	Surface tension for two phase water/steam as a function of T

		Vapour fraction
			x_ph	Vapour fraction as a function of pressure and enthalpy
			x_ps	Vapour fraction as a function of pressure and entropy
		
		Vapour volume fraction
			vx_ph	Vapour volume fraction as a function of pressure and enthalpy
			vx_ps	Vapour volume fraction as a function of pressure and entropy