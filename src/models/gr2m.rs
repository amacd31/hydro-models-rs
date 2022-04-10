#[repr(C)]
pub struct GR2MParams {
    pub x1: f64,
    pub x2: f64,
}

#[repr(C)]
pub struct GR2MModel {
    pub params: GR2MParams,
    pub production_store: f64,
    pub routing_store: f64,
}

impl GR2MModel {
    pub fn create(params: GR2MParams) -> GR2MModel {
        GR2MModel {
            params: params,

            // Completely dry initial catchment
            production_store: 0.,
            routing_store: 0.,
        }
    }

    /*
        Generate monthly simulated streamflow for given rainfall and potential evaporation.

        :param precip: Catchment average rainfall.
        :param potential_evap: Catchment average potential evapotranspiration.

        References:
            Mouelhi, S., Michel, C., Perrin, C., & Andréassian, V. (2006).
                Stepwise development of a two-parameter monthly water
                balance model. Journal of Hydrology, 318(1-4), 200-214.
                http://doi.org/10.1016/j.jhydrol.2005.06.014

            Operation of GR2M:
                http://webgr.irstea.fr/modeles/mensuel-gr2m/fonctionnement-gr2m/?lang=en
    */
    pub fn run(&mut self, precip: &[f64], potential_evap: &[f64]) -> Vec<f64> {
        let mut qsim: Vec<f64> = Vec::with_capacity(precip.len());

        let x1 = self.params.x1;
        let x2 = self.params.x2;

        for (p, e) in precip.iter().zip(potential_evap) {
            let phi = (p / x1).tanh();
            let psi = (e / x1).tanh();

            let s1 =
                (self.production_store + x1 * phi) / (1. + phi * (self.production_store / x1));

            let p1 = p + self.production_store - s1;

            let s2 = s1 * (1. - psi) / (1. + psi * (1. - s1 / x1));

            self.production_store = s2 / (1. + (s2 / x1).powf(3.)).powf(1. / 3.);

            let p2 = s2 - self.production_store;

            let p3 = p1 + p2;

            let r1 = self.routing_store + p3;

            let r2 = x2 * r1;

            let q = r2.powf(2.) / (r2 + 60.);
            qsim.push(q);

            self.routing_store = r2 - q;
        }
        qsim
    }
}