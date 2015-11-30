"""
This module contains all of the cntrl namelist variables for the
different amber programs and automatically loads those dictionaries
with the default values found in that program (sander or pmemd).   

                           GPL LICENSE INFO                             

Copyright (C) 2010  Jason Swails and Ben Roberts

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
   
You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330,
Boston, MA 02111-1307, USA.
"""

class cntrl:

   def __init__(self):
      self.sander  = {'irest': 0, 'ibelly' : 0, 'ntx' : 1, 'ntxo' : 1, 
             'ntcx' : 0, 'ig' : 71277, 'tempi' : 0.0, 'ntb' : -1, 'ntt' : 0, 
             'nchain' : 1, 'temp0' : 300.0, 'tautp' : 1.0, 'ntp' : 0, 
             'pres0' : 1.0, 'comp' : 44.6, 'taup' : 1.0, 'nscm' : 1000, 
             'nstlim' : 1, 'dt' : 0.001, 'ntc' : 1, 'ntcc' : 0, 'nconp' : 0, 
             'tol' : 0.00001, 'ntf' : 1, 'nsnb' : 25, 'cut' : 0.0, 'dielc' : 0,
             'ntpr' : 50, 'ntwx' : 0, 'ntwv' : 0, 'ntwe' : 0, 'ntave' : 0,
             'ntpp' : 0, 'ioutfm' : 0, 'ntr' : 0, 'nrc' : 0, 'ntrx' : 1,
             'taur' : 0, 'nmropt' : 0, 'ivcap' : 0, 'cutcap' : 0.0, 
             'xcap' : 0.0, 'ycap' : 0.0, 'zcap' : 0.0, 'fcap' : 1.5, 
             'xlorth' : -1.0, 'ylorth' : -1.0, 'zlorth' : -1.0, 
             'xorth' : 47114711.0, 'yorth' : 47114711.0, 'zorth' : 47114711.0,
             'forth' : 1.5, 'imin' : 0, 'drms' : 1.0e-4, 'dele' : 0,
             'dx0' : 0.01, 'pencut' : 0.1, 'ipnlty' : 1, 'iscale' : 0,
             'scalm' : 100.0, 'noeskp' : 1, 'maxcyc' : 1, 'ncyc' : 10, 
             'ntmin' : 1, 'vlimit' : 20.0, 'mxsub' : 1, 'ipol' : 0, 
             'jfastw' : 0, 'watnam' : '    ', 'owtnm' : '    ', 
             'hwtnm1' : '    ', 'hwtnm2' : '    ', 'iesp' : 0, 'skmin' : 50, 
             'skmax' : 100, 'vv' : 0, 'vfac' : 0, 'tmode' : 1, 'ips' : 0, 
             'mipsx' : -1, 'mipsy' : -1, 'mipsz' : -1, 'mipso' : 4, 
             'gridips' : 2, 'raips' : -1.0, 'dvbips' : 1.0e-8, 'isgld' : 0, 
             'isgsta' : 1, 'isgend' : 0, 'tsgavg' : 0.2, 'sgft' : 0.0, 
             'tempsg' : 1.0, 'jar' : 0, 'iamoeba' : 0, 'numexchg' : 0, 
             'repcrd' : 1, 'numwatkeep' : -1, 'hybridgb' : 0, 'ntwprt' : 0, 
             'tausw' : 0.1, 'ntwr' : 500, 'iyammp' : 0, 'imcdo' : -1,
             'igb' : 0, 'alpb' : 0, 'arad' : 15.0, 'rgbmax' : 25.0, 
             'saltcon' : 0.0, 'offset' : -999999.0, 'ntwf' : 0,
             'gbsa' : 0, 'vrand' : 1000, 'surften' : 0.005, 'iwrap' : 0, 
             'nrespa' : 1, 'nrespai' : 1, 'gamma_ln' : 0.0, 'extdiel' : 78.5,
             'intdiel' : 1.0, 'cut_inner' : 8.0, 'icfe' : 0, 'clambda' : 0.0,
             'klambda' : 1, 'rbornstat' : 0, 'lastrst' : 1, 'lastist' : 1,
             'itgtmd' : 0, 'tgtrmsd' : 0, 'tgtmdfrc' : 0, 'tgtfitmask' : '',
             'tgtrmsmask' : '', 'idecomp' : 0, 'temp0les' : -1.0, 
             'restraintmask' : '', 'restraint_wt' : 0.0, 'bellymask' : '', 
             'noshakemask' : '', 'crgmask' : '', 'iwrap_mask' : '',
             'mmtsb_switch' : 0, 'mmtsb_iterations' : 100, 'rdt' : 0.0, 
             'icnstph' : 0, 'solvph' : 7.0, 'ntcnstph' : 10, 'ifqnt' : 0, 
             'ievb' : 0, 'ipimd' : 0, 'itimass' : 0, 'ineb' : 0,
             'profile_mpi' : 0, 'ilscivr' : 0, 'icorf_lsc' : 0, 'ipb' : 0, 
             'inp' : 2, 'gbneckscale' : -999999.0, 'gbalphah' : 0.788440, 
             'gbbetah' : 0.798699, 'gbgammah' : 0.437334,  
             'gbalphac' : 0.733756, 'gbbetac' : 0.506378, 
             'gbgammac' : 0.205844, 'gbalphan' : 0.503364, 
             'gbbetan' : 0.316828, 'gbgamman' : 0.192915, 
             'gbalphaos' : 0.867814, 'gbbetaos' : 0.876635,
             'gbgammaos' : 0.387882, 'gbalphap' : 1.0, 'gbbetap' : 0.8, 
             'gbgammap' : 4.851, 'sh' : 1.425952, 'sc' : 1.058554,
             'sn' : 0.733599, 'so' : 1.061039, 'ss' : -0.703469, 'sp' : 0.5,
             't' : 0.0, 'ntn' : 0 , 'scalpha' : 0.5, 'scbeta' : 12.0, 
             'ifsc' : 0, 'scmask' : '', 'logdvdl' : 0, 'dvdl_norest' : 0, 
             'dynlmb' : 0, 'ifmbar' : 0, 'bar_intervall' : 100, 
             'bar_l_min' : 0.1, 'bar_l_max' : 0.9, 'bar_l_incr' : 0.1,
             'idssp' : 0, 'irism' : 0, 'restart_cmd' : '.false.', 
             'eq_cmd' : '.false.', 'adiab_param' : 1.0, 'dec_verbose' : 3,
             'mccycles' : 1, 'ifcr' : 0, 'iamd' : 0, 'barostat' : 1,
             'mcbarint' : 100, 'lj1264' : 0}

      self.pmemd  = {'imin' : 0, 'nmropt' : 0, 'ntx' : 1, 'irest' : 0, 
             'ntrx' : 1, 'ntxo' : 1, 'ntpr' : 50, 'ntave' : 0, 'ntwr' : 500, 
             'iwrap' : 0, 'ntwx' : 0, 'ntwv' : 0, 'ntwe' : 0, 'ioutfm' : 0, 
             'ntwprt' : 0, 'ntf' : 1, 'ntb' : -1, 'dielc' : 1.0, 'cut' : 0.0,
             'nsnb' : 25, 'ipol' : 0, 'igb' : 0, 'intdiel' : 1.0,
             'extdiel' : 78.5, 'saltcon' : 0.0, 'rgbmax' : 25.0,
             'rbornstat' : 0, 'offset' : 0.09, 'gbsa' : 0, 'surften' : 0.005,
             'nrespai' : 1, 'cut_inner' : 8.0, 'ibelly' : 0, 'ntr' : 0,
             'bellymask' : 0, 'maxcyc' : 1, 'ncyc' : 10, 'ntmin' : 1,
             'dx0' : 0.01, 'drms' : 0.0001, 'nstlim' : 1, 'nscm' : 1000,
             't' : 0.0, 'dt' : 0.001, 'nrespa' : 1, 'ntt' : 0, 'temp0' : 300.0,
             'tempi' : 0.0, 'ig' : 71277, 'tautp' : 1.0, 'gamma_ln' : 0.0,
             'vrand' : 1000, 'vlimit' : 20.0, 'ntp' : 0, 'pres0' : 1.0,
             'comp' : 44.6, 'taup' : 1.0, 'ntc' : 1, 'tol' : 0.00001,
             'jfastw' : 0, 'hwtnm1' : '    ', 'hwtnm2' : '    ',
             'owtnm' : '    ', 'watnam' : '    ', 'ivcap' : 0, 'fcap' : 1.5,
             'pencut' : 0.1, 'ndfmin' : 0, 'jar' : 0,
             'no_intermolecular_bonds' : 1, 'ene_avg_sampling' : -1,
             'mdinfo_flush_interval' : 60, 'mdout_flush_interval' : 300,
             'dbg_atom_redistribution' : 0, 'loadbal_verbose' : 0,
             'es_cutoff' : 0.0, 'vdw_cutoff' : 0.0, 'ntwf' : 0,
             'dtemp' : 0, 'dxm' : 0, 'heat' : 0, 'alpb' : 0, 'arad' : 15.0,
             'iamd' : 0, 'restraintwt' : 0.0, 'barostat' : 1, 'mcbarint' : 100,
             'lj1264' : 0, 'icnstph' : 0, 'solvph' : 7.0, 'ntcnstph' : 10}
