import numpy as np
import scipy.constants as const
import pygsl.testing.sf as gsl_sf

from concurrent.futures import ProcessPoolExecutor
from functools import partial,lru_cache

hbar = (const.hbar / (const.elementary_charge*1e6))*(const.c)/(1e-15)
amu = const.physical_constants['atomic mass constant energy equivalent in MeV'][0]
alpha = const.fine_structure

#############################################################################################

def Argonaut(x_points,exc_arr,qvals_dch,gammastates,BR_dch,lstates_dch,m,z_dch,
             sigma_arr,scale_arr,dch_sel=None,r0_user=1.4):
    """
    R-matrix sequential decay strength profile tool. See "Nuclear Reactions for Astrophysics" by Thompson and Nunes for theory overview. 

    Parameters
    ----------

    x_points        : array_like
                        x points over which the R-matrix calculation is being evaluated (in MeV)

    exc_arr         : array_like
                        Array of state excitation energies which will be calculated (in MeV)

    qvals_dch       : array_like
                        An array of the decay channel Q-values (in MeV)

    gammastates     : array_like
                        An array containting the state widths (in MeV)

    BR_dch          : array_like
                        An array of arrays corresponing to the BR of each state into a given decay channel

    lstates_dch     : array_like
                        An array of arrays corresponding to the angular momentum value of each state decaying into the specified decay channel

    m               : array_like
                        Array of arrays corresponding to masses of decay channel constituents (in amu)

    z_dch           : array_like
                        Array of arrays corresponding to the charges of decay channel constituents 

    sigma_arr       : array_like
                        An array of energy dependent standard deviation values used during gaussian convolution (should be same size as x_points)

    scale_arr       : array_like
                        Scale parameters for each state

    dch_sel         : int
                        Selection varaible for only returning voigt lineshape in given decay channel

    r0_user         : float
                        R-matrix channel radius (in fm)

    

    Outputs
    -------
    array_like
        Total R-matrix level profile
    """ 

    @lru_cache(maxsize=1024)
    def _get_coulomb_waves(eta, z, l):
        # cached gsl coloumb wave functions
        return gsl_sf.coulomb_wave_FG_e(eta, z, l, 0)


    def SE(l,eta,z):
        """
        Shift function calculator

        Parameters
        ----------

        l       : float
                    Angular momentum
        eta     : float
                    Sommerfeld parameter
        z       : float
                    Dimensionless radius (z = k*a)

        Outputs
        -------
        float
            Evaluated shift function
        """
        _,F,Fdot,G,_,Gdot,_,_,_,_,_ = _get_coulomb_waves(eta,z,l)   #obtains the F,Fdot,G,Gdot coulomb wavefunctions from cache
        shift = z*(Fdot*F + Gdot*G)/(F**2 + G**2)   # calculates the shift function
        return shift
    

    def rmatrixcalc(engarr,ex,qval,mu_sys,a_sys,zp1,zp2,l_states,gamma_states):
        """
        R-matrix calculation function

        Parameters
        ----------

        engarr      : array_like
                        Array of energy over which the R-matrix calculation needs to be evaluated
        ex          : array_like
                        List of excitation energies of states
        qval        : float
                        Q-value of decay channel (in MeV)
        mu_sys      : float
                        Reduced mass of decay channel constituents
        a_sys       : float
                        R-matrix channel radius (in fm)
        zp1         : float
                        Charge of particle 1 of decay channel constituents
        zp2         : float
                        Charge of particle 2 of decay channel constituents
        l_states    : array_like
                        Angular momentum values of each state in decay
        gamma_states: array_like
                        Widths of each state

        Outputs
        -------
        Array_like
            Evaluated R-matrix level profile
        
        """
        widthtot=[] # init empty lists for storing width and shift 
        shifttot = []
        for n in range(len(ex)):    # loop over each state in excitation array
            if ex[n]+qval > 0:  # if the state is above the decay threshold evaluate normally
                k_resonance = np.sqrt(2*mu_sys*amu*(ex[n]+qval))/hbar   # calculate the wavenumber
                eta_resonance = alpha*zp1*zp2*np.sqrt((mu_sys*amu)/(2*(ex[n]+qval)))    # calculate the sommerfeld parameter
                z_resonance = k_resonance*a_sys # calculate the dimensionless radius

                _,F_res,Fdot_res,G_res,_,Gdot_res,_,_,_,_,_ = _get_coulomb_waves(eta_resonance,z_resonance,l_states[n]) # obtain the Coloumb wavefunctions
                Pl_res = z_resonance/(F_res**2 + G_res**2)  # calculate the penetrability factor
                shiftbc = Pl_res*(Fdot_res*F_res + Gdot_res*G_res)  # calculate the shift funcion boundary condition
                redgamma2_full = gamma_states[n]/(2*Pl_res) # reduced width at pole - ignoring shift correction because of construction of boundary condition (i.e. shift = 0 at resonance)
            else:
                # print(f"Below threshold: {ex[n]+qval}")
                widthtot.append(np.zeros_like(engarr))  # if below the decay threshold set width and shit array to zeros
                shifttot.append(np.zeros_like(engarr))
                continue   
            width_full = [] # init more empty arrays
            shift_full = []
            for eng in engarr:  # loop over each point in the energy array 
                eng+=qval   # add the qvalue to the point
                if eng < 0: 
                    width_full.append(0)    # if below the decay threshold append zero
                    k_full_neg = np.sqrt(2*mu_sys*amu*-eng)/hbar
                    eta_full_neg = alpha*zp1*zp2*np.sqrt((mu_sys*amu)/(2*-eng))
                    z_full_neg = k_full_neg*a_sys
                    shift_full.append(redgamma2_full*(SE(l_states[n],-eta_full_neg,z_full_neg) - shiftbc))  # calculates the shift function at points below the decay threshold
                    continue
                else:
                    k_full = np.sqrt(2*mu_sys*amu*eng)/hbar #wavenumber calculation p=hbar*k
                    eta_full = alpha*zp1*zp2*np.sqrt((mu_sys*amu)/(2*eng)) #sommerfeld coefficient
                    z_full = k_full*a_sys   # dimnesionless radius 

                    _, F, Fdot, G, _, Gdot, _, _, _, _, _ = _get_coulomb_waves(eta_full,z_full,l_states[n])
                    Pl = z_full/(F**2 + G**2)   # calculate the energy dependent penetrability 
                    shift_E = redgamma2_full*(Pl*(Fdot*F + Gdot*G)-shiftbc) # calculate the full shift function
                    width_full.append((2*Pl*redgamma2_full)) #Gamma=2*Pl(E)*gamma2  
                    shift_full.append(shift_E)  # append the shift
                    ## make new width -- entrane channel -- total channel on denominator 
            widthtot.append(width_full) # append each state array into a parent array
            shifttot.append(shift_full)
        return widthtot,shifttot    # return the collected array of arrays
    
    ####

    r0 = r0_user    # important!
    
    ################################

    widtharr = []
    shiftarr = []

    for i in range(len(BR_dch)):    # loop over each decay channel
        
        mu_dch = (m[i][0]*m[i][1])/(m[i][0]+m[i][1])    # calculate the reduced mass
        a_dch = r0*(m[i][0]**(1/3) + m[i][1]**(1/3))    # calculate the channel radius

        gammastates_dch = np.asarray(gammastates)*np.asarray(BR_dch[i]) # obtain the widths of each state relative with the branching fraction


        widtharr_dch,shiftarr_dch = rmatrixcalc(x_points,exc_arr,qvals_dch[i],mu_dch,a_dch, # obtain the widths and shift function from R-matrix calculations
                                   z_dch[i][0],z_dch[i][1],lstates_dch[i],gammastates_dch)

        widtharr.append(np.asarray(widtharr_dch))   # append to arrays and convert to numpy arrays
        shiftarr.append(np.asarray(shiftarr_dch))


    
    widtharr = np.asarray(widtharr)
    shiftarr = np.asarray(shiftarr)
    combinedvoigtarr=np.zeros_like(widtharr[0]) # create empty array for voigt calculation


# Implemented multithreading.. 

## NEED TO UNDERSTAND WHAT THIS IS REALLY DOING... 
    
    process_partial = partial(process_excitation, widtharr=widtharr, shiftarr=shiftarr, x_points=x_points, exc_arr=exc_arr, sigma_arr=sigma_arr, scale_arr=scale_arr) 
    # what this does is create a new function technically which is process_excitation but with only n as the free variable.. 
    # this is so that it can then be fed into the multithreading correctly... 

    with ProcessPoolExecutor() as executor: #multithreading step -- setting up the parallelisation 
        results = list(executor.map(process_partial, range(len(exc_arr))))  #range here is to tell process_excitation how many states it is iterating over -- see above
        # map does what it says and maps the variables to executor in the correct way
        # list is there to store the ProcessPoolExecutor outputs 

    if dch_sel == None:
        combinedvoigtarr = np.sum(results, axis=0)  # flattens array of combined voigts
    else:
        combinedvoigtarr = results[dch_sel]

    return combinedvoigtarr 


def convolution(spectrum, sigma, E, batch_size=1000):
    """Batched convolution with capability of handling of variable kernel size

    Parameters
    ----------

    spectrum        : array_like
                        R-matrix lineshape

    sigma           : array_like
                        Energy dependent gaussian standard deviation 

    E               : array_like
                        Array of energies at which the R-matrix lineshape corresponds

    batch_size      : int
                        Size of batches evaluated in convolution

    Outputs
    -------
    Array_like
        Returned voigt level profiles for each state
    
    """
    n = len(E)  # total size of kernel points
    result = np.zeros_like(spectrum)    # empty numpy array to build from   
    
    for i in range(0, n, batch_size):   # loop over the total kernel points in groups of the batch size
        end_idx = min(i + batch_size, n)    #return the smallest between the two arguments
        batch_slice = slice(i, end_idx) # slice the batch based on this
        
        distances = E[batch_slice, np.newaxis] - E # calculate distances matrix for this batch
        weights = np.exp(-0.5 * (distances / sigma[batch_slice, np.newaxis])**2)# calculate Gaussian for convolution with BW
        
        norm = np.sqrt(2*np.pi) * sigma[batch_slice, np.newaxis]    # normalise Gaussians
        weights /= norm
        
        row_sums = np.trapz(weights, E, axis=1)# ensure normalization is correct for each point
        weights /= row_sums[:, np.newaxis]
        
        result[batch_slice] = np.dot(weights, spectrum)# apply convolution
    
    return result




# for MT to work this has to be outside so that it is pickleable (not entirely sure what this means... )

def process_excitation(n, widtharr, shiftarr, x_points, exc_arr, sigma_arr, scale_arr):
        widthtemp = np.zeros(len(widtharr[0][n]))
        shifttemp = np.zeros(len(shiftarr[0][n]))

        for dch in range(len(widtharr)):
            widthtemp+=widtharr[dch][n]
            shifttemp+=shiftarr[dch][n] # fills empty arrays with relevant width and shift values


            # need to pull out individual width temp rather than total width temp for numerator!

        denom = (x_points - exc_arr[n] + shifttemp)**2 + 0.25 * widthtemp**2    # calculates the denominator of the R-matrix BW 
        # for i,val in enumerate(denom):
        #     if val == 0.0:
        #         print(i)

        lineshape = np.zeros(len(widtharr[0][n]))   # makes an empty array to store final lineshape
        
        for dch in range(len(widtharr)):    # loops over each decay channel
            # need to divide widthtemp by widtharr[dch][n] but for any zero values in widtharr[dch][n] need to just set the answer to zero
            lineshape += widtharr[dch][n] / denom   # calculates the BW for R-matrix

        '''lineshape = np.zeros(len(widtharr[0][n]))

        for dch in range(len(widtharr)):
            lineshape += widtharr[dch][n]**2 / denom'''

        #lineshape = widthtemp**2 / denom

           ## should not be width temp should be individual width... 

        voigtlineshape = convolution(lineshape, sigma_arr, x_points)    # calculate the voigt

        return voigtlineshape * scale_arr[n]    # return the voigt scaled by a specified amount
