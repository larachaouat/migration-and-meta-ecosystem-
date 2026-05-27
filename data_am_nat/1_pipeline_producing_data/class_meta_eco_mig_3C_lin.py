
# This script
# -  model of a 3-compartment ecosystem with migration beteween two ecosystems

# Import necessary libraries
from math import *
import numpy as np
from scipy.integrate import solve_ivp



class two_eco_s_3C_lin:
    """
    two_eco_s_3C_lin: Coupled two-ecosystem seasonal consumer-producer model (linear producers)
    Description
    -----------
    This class implements a mechanistic, time-varying model of two ecosystems (1 and 2),
    each with nutrients (N), primary producers (P) and resident consumers (C1, C2), plus
    a mobile consumer (migrant) that switches between ecosystems. Seasonal forcing is
    applied to producer uptake rates via a piecewise quadratic (two-part) periodic function.
    Consumer feeding follows a Holling type-II functional response. The model supports
    a migration fraction (Fmig) that divides each seasonal period between ecosystems
    and a mode without migration that integrates a single state equation set for a
    given starting ecosystem.
    Primary model features:
    - State variables per ecosystem: nutrient (N), producer biomass (P), resident consumer (C_res).
    - A single mobile consumer whose presence alternates between ecosystems (consumer in eco1 or eco2).
    - Seasonal modulation of producer uptake rates (a1, a2) via sin_seasonality() with configurable
        amplitude, season length, and optional mean-preserving or extremum-preserving coefficient computation.
    - Holling type-II feeding: feeding_rate = F * P / (1 + F * h * P).
    - Recycling of consumer and producer mortality into nutrients via rc and rp.

    State vector ordering
    ---------------------
    All ODE functions and API methods use an 8-component state vector with the following ordering:
            [N1, P1, C1, C, N2, P2, C2, Ca]
    where
    - N1: nutrient in ecosystem 1
    - P1: primary producer biomass in ecosystem 1
    - C1: resident consumer in ecosystem 1 (consumer 1)
    - C: mobile consumer when currently occupying ecosystem 1 (or placeholder; see methods)
    - N2, P2, C2: analogous quantities for ecosystem 2
    - Ca: mobile consumer when currently occupying ecosystem 2
    API / Public methods
    --------------------
    __init__(F, F1, F2, h1, h2, rp, mp, I1, a1, amp1, e1, rc,
                     mc1, mc2, s, f, mc0=None, per_s1=0.5, per_s2=None,
                     Fmig=0.5, I2=None, a2=None, amp2=None, e2=None, mean_norm=None)
    - Initialize model parameters and precompute seasonal coefficient pairs (a_season1, a_season2).
    - Parameters (brief):
            - F, F1, F2: attack/feeding rates (mobile consumer and resident consumers)
            - h1, h2: handling times for the associated feeding interactions
            - rp: producer mortality recycling fraction (to nutrients)
            - mp: producer mortality rate
            - I1, I2: external nutrient input (I2 defaults to I1)
            - a1, a2: mean uptake coefficients for producers (a2 defaults to a1)
            - amp1, amp2: amplitude multipliers for seasonal forcing (amp2 defaults to 0 if not provided)
            - e1, e2: nutrient loss/export rates (e2 defaults to e1)
            - rc: consumer mortality recycling fraction (to nutrients)
            - mc1, mc2, mc0: consumer mortality rates (mc0 defaults to mc1)
            - s: seasonality synchrony parameter (phase relation between ecosystems)
            - f: seasonal frequency (cycles per unit time; period T = 1/f)
            - per_s1, per_s2: fraction of each period defined as "summer" for each ecosystem (per_s2 defaults to 0.5)
            - Fmig: fraction of a period that the mobile consumer spends in one ecosystem before switching
            - mean_norm: if True use fact_season_mean (same mean), otherwise fact_season_extrem (same extrema)
    - Attributes:
            - model parameters are stored as attributes (a1, a2, amp1, amp2, F, F1, F2, h1, h2, rp, rc, mp, mc1, mc2, mc0, I1, I2, e1, e2, s, f, per_s1, per_s2, Fmig)
            - a_season1, a_season2: coefficient pairs [a, b] used by sin_seasonality()
    fact_season_mean(p)
    - Compute [a, b] coefficient pair that yields a seasonal quadratic piecewise function with a specified
        property (constructed in code for the "same mean" option). Input p is the summer fraction (per_s).
    - Returns: list [a, b].
    fact_season_extrem(p)
    - Compute [a, b] coefficient pair for the "same extrema" variant of the piecewise seasonal shape.
    - Returns: list [a, b].
    sin_seasonality(x, per_s, a_season)
    - Evaluate the piecewise quadratic seasonal factor at time x (continuous time).
    - Inputs:
            - x: time (scalar or numpy array)
            - per_s: fraction of period considered summer
            - a_season: coefficient pair [a, b] obtained from fact_season_* methods
    - Behavior:
            - Period T = 1/f (f is model attribute).
            - The function is periodic with period T and is defined by two quadratics:
                    - summer segment (0 <= t_mod < per_s * T): a * t_mod * (t_mod - per_s * T)
                    - winter segment (per_s * T <= t_mod < T): b * (t_mod - T) * (t_mod - per_s * T)
            - Typically used multiplicatively as (1 + sin_seasonality(...) * amp) to modulate uptake.
    set_eq_1(t, A)
    - ODE right-hand side when the mobile consumer is present in ecosystem 1.
    - Inputs:
            - t: time
            - A: length-8 state vector in the ordering described above
    - Returns:
            - list of 8 derivatives in the same ordering: [dN1dt, dP1dt, dC1dt, dCdt, dN2dt, dP2dt, dC2dt, dCadt]
    - Notes:
            - Feeding terms use Holling type-II forms with appropriate F, F1, F2 and handling times.
            - Seasonality applied to uptake terms for a1 and a2 with user-specified phase offsets via s.
    set_eq_2(t, A)
    - ODE right-hand side when the mobile consumer is present in ecosystem 2.
    - Same signature and returned ordering as set_eq_1, but with the mobile consumer affecting ecosystem 2.
    solve_model(X0, nb_years, eco)
    - Integrate the model with migration (consumer alternates between ecosystems inside each period).
    - Inputs:
            - X0: initial length-8 state vector ordered as [N1, P1, C1, C, N2, P2, C2, Ca]
            - nb_years: number of annual cycles to simulate (integration proceeds in seasonal pieces)
            - eco: integer 1 or 2 indicating which ecosystem the mobile consumer occupies at simulation start
    - Behavior:
            - The method alternates calls to set_eq_1 and set_eq_2 splitting each seasonal period by Fmig
                (the mobile consumer's residence time fraction).
            - Integrates in pieces using scipy.integrate.solve_ivp (max_step is set in code).
            - Concatenates piecewise solutions and returns full trajectories.
    - Returns:
            - time: 1D numpy array of time points
            - N1, P1, C1, C, N2, P2, C2, Ca: 1D numpy arrays each aligned to 'time'
    solve_model_nomig(X0, nb_years, eco)
    - Integrate the model without migration: a single set of ODEs (either set_eq_1 or set_eq_2)
        is integrated continuously for nb_years * (1/f) time units.
    - Inputs and returns similar to solve_model.
    Dependencies
    ------------
    - numpy (as np) for vectorized seasonal evaluation
    - scipy.integrate.solve_ivp for numerical integration

  
    """
    def __init__(self,   F, F1, F2, h1, h2,  rp,  mp, I1, a1,  amp1,  e1, rc,  mc1, mc2, s, f, mc0= None, per_s1=.5, per_s2= None, Fmig=.5, I2= None, a2= None, amp2= None,  e2=None, mean_norm=None ):
        self.a1 = a1  # mean uptake primary producer 1
        self.a2 = a2 if a2 is not None else a1  # mean uptake primary producer 2
        
        self.amp1 = amp1  # amplitude
        self.amp2 = amp2 if amp2 is not None else 0 # amplitude

        self.F = F    #attack rate migrant
        self.F1 = F1  # attack rate
        self.F2 = F2  # attack rate
        self.h1 = h1  # handling time
        self.h2 = h2  # handling time

        self.rp = rp  # recycle producer
        self.rc = rc  # recycle consumer


        self.I1 = I1  # input
        self.I2 = I2 if I2 is not None else I1  # input

        self.e1 = e1  # output
        self.e2 = e2 if e2 is not None else e1


        self.mp = mp #mortality producer
        self.mc1 = mc1 #mortality consumer
        self.mc2 = mc2 #mortality consumer
        self.mc0 = mc0 if mc0 is not None else mc1 #mortality conso in 2

        self.s = s #synchronicity of seasons
        self.f = f #frequency season
 
        self.Fmig = Fmig #migration threshold


        self.per_s1 = per_s1 #percentage of summer in one year
        self.per_s2 = per_s2 if per_s2 is not None else .5 #percentage of summer in one year
        self.a_season1 = self.fact_season_mean( self.per_s1) if mean_norm== True else self.fact_season_extrem( self.per_s1)
        self.a_season2 = self.fact_season_mean( self.per_s2) if mean_norm== True else self.fact_season_extrem( self.per_s2)

        self.Fmig = Fmig  # time spent in eco 1 vs time spent in eco 2

    ## same mean
    def fact_season_mean(self,p):
        if p >= .5:
            alpha = (p/self.f+1/self.f)/2
            b= -1/((alpha-p/self.f)*(alpha-1/self.f))    
            a = b* ((p/self.f)**3-3*(p**2/self.f**3)+3*(p/self.f**3)-(1/self.f)**3)/(p/self.f)**3

        else: 
            alpha = (p/self.f+0)/2
            a= 1/(alpha*(alpha-p/self.f))
            b = a/((p/self.f)**3-3*(p**2/self.f**3)+3*(p/self.f**3)-(1/self.f)**3)*(p/self.f)**3

        return [a,b]
    ## same extrema
    def fact_season_extrem(self,p):
        alpha = (p/self.f+1/self.f)/2
        b= -1/((alpha-p/self.f)*(alpha-1/self.f))    
        alpha = (p/self.f+0)/2
        a = 1/((alpha-p/self.f)*(alpha))  

        return [a,b]
    def sin_seasonality(self, x, per_s, a_season ):
        T=1/self.f
        x_normalized = (x+self.per_s1*T) % T
        a= a_season[0]
        b= a_season[1]

        y = np.where(x_normalized < T*per_s,
                    a*(x_normalized)*(x_normalized-per_s*T),  # Positive part  = summer
                    b*(x_normalized-T)*(x_normalized-per_s*T)) # Negative part  = winter
        return y
     #consumers in ecosystem 1
    def set_eq_1(self, t, A):
        N1, P1, C1, C, N2, P2, C2, Ca = A
        dCdt = -self.mc1*C + (self.F * P1 / (1 + self.F * self.h1  * P1 )) * C  
        dCadt = 0

        dN1dt = self.I1- self.e1 * N1 -(self.a1*(1+self.sin_seasonality(t+1/(self.f*self.s), self.per_s1, self.a_season1)*self.amp1)) * N1 * P1 + self.rc*(self.mc1 * C + self.mc2 * C1) + self.rp*self.mp * (P1)
        dP1dt = (self.a1*(1+self.sin_seasonality(t+1/(self.f*self.s), self.per_s1, self.a_season1)*self.amp1)) * N1 * P1  - self.mp*P1 - self.F  * P1 / (1 + self.F * self.h1 * P1 ) * C - (self.F1 * P1 / (1 + self.F1 * self.h2  * P1 )) * C1
        dC1dt = -self.mc2*C1 + (self.F1 * P1 / (1 + self.F1 * self.h2  * P1 )) * C1

        dN2dt = self.I2- self.e2 * N2 -(self.a2*(1+self.sin_seasonality(t+1/(self.f), self.per_s2, self.a_season2)*self.amp2)) * N2 * P2 + self.rp*self.mp * (P2)+ self.rc*self.mc2 * C2
        dP2dt = (self.a2*(1+self.sin_seasonality(t, self.per_s2, self.a_season2)*self.amp2)) * N2  * P2  - self.mp*P2 - self.F2 * P2 / (1 + self.F2 * self.h2 *  P2 ) * ( C2)
        dC2dt = -self.mc2*C2 + (self.F2 * P2 / (1 + self.F2 * self.h2  * P2 )) * C2
        return [dN1dt, dP1dt, dC1dt, dCdt, dN2dt,dP2dt, dC2dt, dCadt]
    
    #consumers in ecosystem 2
    def set_eq_2(self, t, A):
        N1, P1, C1, C, N2, P2, C2, Ca= A
        dCdt = 0
        dCadt = -self.mc0*Ca + (self.F * P2 / (1 + self.F * self.h1  * P2 )) * Ca 
        
        dN1dt = self.I1 - self.e1 * N1 - (self.a1*(1+self.sin_seasonality(t+1/(self.f*self.s), self.per_s1, self.a_season1)*self.amp1)) * N1  * P1 + self.rc*(self.mc2 * C1) + self.rp*self.mp * (P1)
        dP1dt = (self.a1*(1+self.sin_seasonality(t+1/(self.f*self.s), self.per_s1, self.a_season1)*self.amp1)) * N1 * P1  - self.mp*P1 -(self.F1 * P1 / (1 + self.F1 * self.h2  * P1 )) * C1
        dC1dt = -self.mc2*C1 + (self.F1 * P1 / (1 + self.F1 * self.h2  * P1 )) * C1

        dN2dt = self.I2 - self.e2 * N2 - (self.a2*(1+self.sin_seasonality(t, self.per_s2, self.a_season2)*self.amp2)) * N2 * P2 + self.rp*self.mp * (P2) + self.rc*(self.mc0 * Ca + self.mc2 * C2)
        dP2dt = (self.a2*(1+self.sin_seasonality(t, self.per_s2, self.a_season2)*self.amp2)) * N2 * P2  - self.mp*P2 - self.F * P2 / (1 + self.F * self.h1 *  P2 ) * Ca  -self.F2 * P2 / (1 + self.F2 * self.h2 *  P2 ) * C2

        dC2dt = -self.mc2*C2 + (self.F2 * P2 / (1 + self.F2 * self.h2  * P2 )) * C2
        return [dN1dt, dP1dt, dC1dt, dCdt, dN2dt,dP2dt, dC2dt, dCadt]


    def solve_model(self, X0, nb_years, eco):
        state = [self.set_eq_1, self.set_eq_2]
        ts = []
        ys = []
        t = 0
        X=X0
        state_i= eco-1 #eco2
        for i in range(1, 2*(nb_years+1)): # 2 seasons per year
                if state_i == 0: #eco1 to 2
                        sol=solve_ivp(state[state_i], [t,t+1/self.f*self.Fmig], [X[0],X[1],X[2], X[7], X[4],X[5],X[6], X[3]], max_step=.1) 
                        state_i =1
                else :
                        sol=solve_ivp(state[state_i], [t,t+1/self.f*(1-self.Fmig)],[X[0],X[1],X[2], X[7], X[4],X[5],X[6], X[3]], max_step=.1) 
                        state_i =0
                if t < 1000:
                    ys.append(sol.y[:,:-1])
                    ts.append(sol.t[:-1])
                else:
                    ys.append(sol.y[:,:-2])
                    ts.append(sol.t[:-2])
                # New start time for integration
                t = sol.t[-1]
                # Reset initial state
                X = sol.y[:, -1].copy()
        time = np.concatenate(ts)
        y = np.concatenate(ys, axis=1).T
        N1, P1, C1, C, N2, P2, C2 , Ca = y[:, :8].T
        
        return time, N1, P1, C1, C, N2, P2, C2, Ca
    
    def solve_model_nomig (self, X0, nb_years, eco):
        state = [self.set_eq_1, self.set_eq_2]
        X=X0
        state_i= eco-1
        sol=solve_ivp(state[state_i], [0,1/self.f*nb_years], [X[0],X[1],X[2], X[7], X[4],X[5],X[6], X[3]], max_step=.1) 

        time = sol.t
        N1, P1, C1, C, N2, P2, C2, Ca = sol.y
        
        return time, N1, P1, C1, C, N2, P2, C2, Ca


