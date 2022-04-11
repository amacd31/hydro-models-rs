use std::collections::HashMap;

#[cfg(feature = "python")]
use pyo3::prelude::*;


/*
Unit hydrograph ordinates for UH1 derived from S-curves.
*/
pub(crate) fn s_curves1(t: f64, x4: f64) -> f64 {
    match t {
        t if t <= 0. => 0.,
        t if t < x4 => (t / x4).powf(2.5),
        // t >= x4
        _ => 1.,
    }
}

/*
Unit hydrograph ordinates for UH2 derived from S-curves.
*/
pub(crate) fn s_curves2(t: f64, x4: f64) -> f64 {
    match t {
        t if t <= 0. => 0.,
        t if t < x4 => 0.5 * (t / x4).powf(2.5),
        t if t < 2. * x4 => 1. - 0.5 * (2. - t / x4).powf(2.5),
        // t >= 2*x4
        _ => 1.,
    }
}

#[cfg_attr(feature = "python", pyclass)]
#[derive(Clone)]
#[repr(C)]
pub struct GR4JParams {
    pub x1: f64,
    pub x2: f64,
    pub x3: f64,
    pub x4: f64,
}

#[cfg_attr(feature = "python", pyclass)]
#[repr(C)]
pub struct GR4JModel {
    pub params: GR4JParams,
    pub production_store: f64,
    pub routing_store: f64,
    pub uh1: Vec<f64>,
    pub uh2: Vec<f64>,
}

impl Default for GR4JModel {
    fn default() -> GR4JModel {
        GR4JModel {
            // Default params are the median values from
            // Perrin, Charles, Claude Michel, and Vazken Andréassian.
            // "Improvement of a parsimonious model for streamflow simulation."
            // Journal of Hydrology 279, no. 1 (2003): 275-289.
            params: GR4JParams {
                x1: 350.,
                x2: 0.,
                x3: 90.,
                x4: 1.7,
            },

            // Completely dry initial catchment
            production_store: 0.,
            routing_store: 0.,

            // Unit hydrographs initialized to zero and size set based on x4
            uh1: vec![0.; (1.7f64).ceil() as usize],
            uh2: vec![0.; (2.0 * 1.7f64).ceil() as usize],
        }
    }
}

impl GR4JModel {
    pub fn create(params: GR4JParams) -> GR4JModel {
        let n_uh1 = params.x4.ceil() as usize;
        let n_uh2 = (2.0 * params.x4).ceil() as usize;

        // Default unit hydrographs to all zeroes when initializing
        let uh1 = vec![0.; n_uh1];
        let uh2 = vec![0.; n_uh2];

        GR4JModel {
            params: params,

            // Completely dry initial catchment
            production_store: 0.,
            routing_store: 0.,

            uh1: uh1,
            uh2: uh2,
        }
    }

    /*
        Generate simulated streamflow for given rainfall and potential evaporation.

        The resulting simulation is appended to the vector stored in the qsim field.

        :param precip: Catchment average rainfall.
        :param potential_evap: Catchment average potential evapotranspiration.
    */
    pub fn run(&mut self, precip: &[f64], potential_evap: &[f64]) -> Vec<f64> {
        let mut qsim: Vec<f64> = Vec::with_capacity(precip.len());

        let n_uh1 = self.params.x4.ceil() as i32;
        let n_uh2 = (2.0 * self.params.x4).ceil() as i32;

        let mut uh1_ordinates = vec![0.; n_uh1 as usize];
        let mut uh2_ordinates = vec![0.; n_uh2 as usize];

        for t in 1..(n_uh1 + 1) {
            uh1_ordinates[(t - 1) as usize] = s_curves1(f64::from(t), self.params.x4)
                - s_curves1(f64::from(t) - 1., self.params.x4);
        }

        for t in 1..(n_uh2 + 1) {
            uh2_ordinates[(t - 1) as usize] = s_curves2(f64::from(t), self.params.x4)
                - s_curves2(f64::from(t) - 1., self.params.x4);
        }

        for (p, e) in precip.iter().zip(potential_evap) {
            let net_evap;
            let mut routing_pattern;
            let reservoir_production;
            if p > e {
                net_evap = 0.;
                let scaled_net_precip = (13.0f64).min((p - e) / self.params.x1);
                let tanh_scaled_net_precip = scaled_net_precip.tanh();
                reservoir_production = (self.params.x1
                    * (1. - (self.production_store / self.params.x1).powf(2.))
                    * tanh_scaled_net_precip)
                    / (1. + self.production_store / self.params.x1 * tanh_scaled_net_precip);

                routing_pattern = p - e - reservoir_production;
            } else {
                let scaled_net_evap = (13.0f64).min((e - p) / self.params.x1);
                let tanh_scaled_net_evap = scaled_net_evap.tanh();

                let ps_div_x1 =
                    (2. - self.production_store / self.params.x1) * tanh_scaled_net_evap;
                net_evap = self.production_store * (ps_div_x1)
                    / (1.
                        + (1. - self.production_store / self.params.x1) * tanh_scaled_net_evap);

                reservoir_production = 0.;
                routing_pattern = 0.;
            }

            self.production_store = self.production_store - net_evap + reservoir_production;

            let percolation = self.production_store
                / (1. + (self.production_store / 2.25 / self.params.x1).powf(4.)).powf(0.25);
            routing_pattern += self.production_store - percolation;

            self.production_store = percolation;

            for i in 0..(self.uh1.len() - 1) {
                self.uh1[i] = self.uh1[i + 1] + uh1_ordinates[i] * routing_pattern;
            }
            if let (Some(last_uh1), Some(last_ordinate)) =
                (self.uh1.last_mut(), uh1_ordinates.last())
            {
                *last_uh1 = *last_ordinate * routing_pattern;
            }

            for j in 0..(self.uh2.len() - 1) {
                self.uh2[j] = self.uh2[j + 1] + uh2_ordinates[j] * routing_pattern
            }
            if let (Some(last_uh2), Some(last_ordinate)) =
                (self.uh2.last_mut(), uh2_ordinates.last())
            {
                *last_uh2 = *last_ordinate * routing_pattern;
            }

            let groundwater_exchange =
                self.params.x2 * (self.routing_store / self.params.x3).powf(3.5);
            self.routing_store =
                (0.0f64).max(self.routing_store + self.uh1[0] * 0.9 + groundwater_exchange);

            let r2 = self.routing_store
                / (1. + (self.routing_store / self.params.x3).powf(4.)).powf(0.25);
            let qr = self.routing_store - r2;
            self.routing_store = r2;
            let qd = (0.0f64).max(self.uh2[0] * 0.1 + groundwater_exchange);
            let q = qr + qd;

            qsim.push(q);
        }
        qsim
    }

    pub fn init(
        &mut self,
        params: &HashMap<&str, f64>,
        production_store: Option<f64>,
        routing_store: Option<f64>,
        unit_hydrographs: Option<HashMap<&str, Vec<f64>>>,
    ) {
        self.params.x1 = params["X1"];
        self.params.x2 = params["X2"];
        self.params.x3 = params["X3"];
        self.params.x4 = params["X4"];

        if let Some(ps) = production_store {
            self.production_store = ps;
        }

        if let Some(rs) = routing_store {
            self.routing_store = rs;
        }

        let n_uh1 = self.params.x4.ceil() as usize;
        let n_uh2 = (2.0 * self.params.x4).ceil() as usize;

        // Default to all zeroes but set after if defined in unit_hydrographs HashMap
        self.uh1 = vec![0.; n_uh1];
        self.uh2 = vec![0.; n_uh2];

        if let Some(unit_hydrographs) = unit_hydrographs {
            if let Some(uh1) = unit_hydrographs.get("uh1") {
                self.uh1.clone_from(&*uh1);
            }

            if let Some(uh2) = unit_hydrographs.get("uh2") {
                self.uh2.clone_from(&*uh2);
            }
        }
    }
}
