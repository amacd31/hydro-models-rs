pub mod hydromodels {
    use std::collections::HashMap;

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

    #[repr(C)]
    pub struct GR4JParams {
        pub x1: f64,
        pub x2: f64,
        pub x3: f64,
        pub x4: f64,
    }

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
}

#[cfg(test)]
mod tests {
    use hydromodels;
    use std::collections::HashMap;

    #[test]
    fn gr4j_create_test() {
        let params = hydromodels::GR4JParams {
            x1: 10.,
            x2: 5.,
            x3: 4.,
            x4: 1.,
        };

        let expected = vec![
            0.13030636843636356,
            0.7348755690383941,
            3.8482364547176102,
            11.246676991863168,
            13.90322269162079,
        ];

        let mut gr4j = hydromodels::GR4JModel::create(params);

        let qsim = gr4j.run(&[10., 2., 3., 4., 5.], &[0.5, 0.5, 0.5, 0.5, 0.5]);

        assert_eq!(qsim, expected);
    }

    #[test]
    fn gr2m_test() {
        let params = hydromodels::GR2MParams { x1: 1., x2: 1. };

        let expected = vec![
            0.0009450053530675536,
            0.0327630181266177,
            0.2055971698811445,
            0.6632149029988136,
            0.991927371887451,
            1.0533436405957102,
            0.9029856114321761,
        ];

        let mut gr2m = hydromodels::GR2MModel::create(params);

        let qsim = gr2m.run(&[1., 2., 3., 4., 3., 2., 1.], &[1., 1., 1., 1., 1., 1., 1.]);

        assert_eq!(gr2m.production_store, 0.1792526700466929);
        assert_eq!(gr2m.routing_store, 6.922989036389427);
        assert_eq!(qsim, expected);
    }

    #[test]
    fn gr4j_simple_test() {
        let mut params = HashMap::new();

        params.insert("X1", 10.);
        params.insert("X2", 5.);
        params.insert("X3", 4.);
        params.insert("X4", 1.);

        let expected = vec![
            0.13030636843636356,
            0.7348755690383941,
            3.8482364547176102,
            11.246676991863168,
            13.90322269162079,
        ];

        let mut gr4j = hydromodels::GR4JModel {
            ..Default::default()
        };

        gr4j.init(&params, None, None, None);
        let qsim = gr4j.run(&[10., 2., 3., 4., 5.], &[0.5, 0.5, 0.5, 0.5, 0.5]);

        assert_eq!(qsim, expected);
    }

    #[test]
    fn gr4j_scaled_test() {
        let mut params = HashMap::new();

        params.insert("X1", 10.);
        params.insert("X2", 5.);
        params.insert("X3", 4.);
        params.insert("X4", 1.);

        let expected = vec![
            0.13030636843636356,
            0.5196067786804051,
            0.7978026532390882,
            128.2406932741725,
            20.992238476264593,
        ];

        let mut gr4j = hydromodels::GR4JModel {
            ..Default::default()
        };

        gr4j.init(&params, None, None, None);
        let qsim = gr4j.run(&[10., 2., 3., 150., 5.], &[0.5, 14., 0.5, 10., 0.5]);

        assert_eq!(qsim, expected);
    }

    #[test]
    fn gr4j_set_production_store() {
        let mut params = HashMap::new();

        params.insert("X1", 10.);
        params.insert("X2", 5.);
        params.insert("X3", 4.);
        params.insert("X4", 1.);

        let mut states = HashMap::new();

        states.insert("production_store", 10.);

        let expected = vec![
            5.1602235324393675,
            10.091188499285725,
            9.82974398339987,
            136.95699908376093,
            21.019904684254975,
        ];

        let mut gr4j = hydromodels::GR4JModel {
            ..Default::default()
        };

        gr4j.init(&params, Some(10.), None, None);
        let qsim = gr4j.run(&[10., 2., 3., 150., 5.], &[0.5, 14., 0.5, 10., 0.5]);

        assert_eq!(qsim, expected);
    }

    #[test]
    fn gr4j_set_unit_hydrographs() {
        let mut params = HashMap::new();

        params.insert("X1", 10.);
        params.insert("X2", 5.);
        params.insert("X3", 4.);
        params.insert("X4", 2.);

        let mut states = HashMap::new();

        states.insert("production_store", 10.);

        let expected = vec![
            7.031781387527497,
            15.758863927257103,
            10.450019503728232,
            32.38841274927161,
            115.58026737087125,
        ];

        let mut unit_hydrographs = HashMap::new();
        unit_hydrographs.insert("uh1", vec![111.13119599074196, 3.7349877368581]);
        unit_hydrographs.insert(
            "uh2",
            vec![
                55.65159496199312,
                57.05053793090891,
                13.713382653511243,
                0.40102046754686754,
            ],
        );

        let mut gr4j = hydromodels::GR4JModel {
            ..Default::default()
        };

        gr4j.init(&params, Some(10.), None, Some(unit_hydrographs));
        let qsim = gr4j.run(&[10., 2., 3., 150., 5.], &[0.5, 14., 0.5, 10., 0.5]);

        assert_eq!(qsim, expected);
    }

    #[test]
    fn gr4j_set_unit_hydrograph_1() {
        let mut params = HashMap::new();

        params.insert("X1", 10.);
        params.insert("X2", 5.);
        params.insert("X3", 4.);
        params.insert("X4", 2.);

        let mut states = HashMap::new();

        states.insert("production_store", 10.);

        let expected = vec![
            1.326727594436605,
            14.387525661905979,
            10.409917456973545,
            32.38841274927161,
            115.58026737087125,
        ];

        let mut unit_hydrographs = HashMap::new();
        unit_hydrographs.insert("uh1", vec![111.13119599074196, 3.7349877368581]);

        let mut gr4j = hydromodels::GR4JModel {
            ..Default::default()
        };

        gr4j.init(&params, Some(10.), None, Some(unit_hydrographs));
        let qsim = gr4j.run(&[10., 2., 3., 150., 5.], &[0.5, 14., 0.5, 10., 0.5]);

        assert_eq!(qsim, expected);
    }

    #[test]
    fn s_curves1_test() {
        assert_eq!(hydromodels::s_curves1(-1.0, 2.0), 0.0);
        assert_eq!(hydromodels::s_curves1(0.0, 2.0), 0.0);
        assert_eq!(hydromodels::s_curves1(1.0, 2.0), 0.1767766952966369);
        assert_eq!(hydromodels::s_curves1(2.0, 2.0), 1.0);
        assert_eq!(hydromodels::s_curves1(2.1, 2.0), 1.0);
    }

    #[test]
    fn s_curves2_test() {
        assert_eq!(hydromodels::s_curves2(-1.0, 2.0), 0.0);
        assert_eq!(hydromodels::s_curves2(0.0, 2.0), 0.0);
        assert_eq!(hydromodels::s_curves2(1.0, 2.0), 0.08838834764831845);
        assert_eq!(hydromodels::s_curves2(2.0, 2.0), 0.5);
        assert_eq!(hydromodels::s_curves2(2.1, 2.0), 0.5601759051904955);
        assert_eq!(hydromodels::s_curves2(4.0, 2.0), 1.0);
        assert_eq!(hydromodels::s_curves2(4.1, 2.0), 1.0);
    }
}
