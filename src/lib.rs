pub mod hydromodels {
    use std::collections::HashMap;

    /*
    Unit hydrograph ordinates for UH1 derived from S-curves.
    */
    pub(crate) fn s_curves1(t: f64, x4: f64) -> f64 {

        if t <= 0.0 {
            return 0.0;
        } else if t < x4 {
            return (t / x4).powf(2.5);
        } else {
            // t >= x4
            return 1.0;
        }
    }

    /*
    Unit hydrograph ordinates for UH2 derived from S-curves.
    */
    pub(crate) fn s_curves2(t: f64, x4: f64) -> f64 {

        if t <= 0.0 {
            return 0.0;
        } else if t < x4 {
            return 0.5 * (t / x4).powf(2.5);
        } else if t < 2.0 * x4 {
            return 1.0 - 0.5 * (2.0 - t / x4).powf(2.5);
        } else {
            // t >= x4
            return 1.0;
        }
    }

    pub struct GR4JModel {
        pub qsim: Vec<f64>,
        pub x1: f64,
        pub x2: f64,
        pub x3: f64,
        pub x4: f64,
        pub production_store: f64,
        pub routing_store: f64,
        pub uh1: Vec<f64>,
        pub uh2: Vec<f64>,
    }

    impl Default for GR4JModel {
        fn default() -> GR4JModel {
            GR4JModel {
                qsim: Vec::new(),
                // Default params are the median values from
                // Perrin, Charles, Claude Michel, and Vazken Andréassian.
                // "Improvement of a parsimonious model for streamflow simulation."
                // Journal of Hydrology 279, no. 1 (2003): 275-289.
                x1: 350.,
                x2: 0.,
                x3: 90.,
                x4: 1.7,

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
        /*
            Generate simulated streamflow for given rainfall and potential evaporation.

            The resulting simulation is appended to the vector stored in the qsim field.

            :param precip: Catchment average rainfall.
            :param potential_evap: Catchment average potential evapotranspiration.
        */
        pub fn run(&mut self, precip: &[f64], potential_evap: &[f64]) {

            let n_uh1 = self.x4.ceil() as i32;
            let n_uh2 = (2.0 * self.x4).ceil() as i32;

            let mut uh1_ordinates: Vec<f64> = Vec::with_capacity(n_uh1 as usize);
            for _i in 0..n_uh1 {
                uh1_ordinates.push(0.);
            }

            let mut uh2_ordinates: Vec<f64> = Vec::with_capacity(n_uh2 as usize);
            for _i in 0..n_uh2 {
                uh2_ordinates.push(0.);
            }

            for t in 1..(n_uh1 + 1) {
                uh1_ordinates[(t - 1) as usize] = s_curves1(t as f64, self.x4) -
                    s_curves1(t as f64 - 1., self.x4);
            }

            for t in 1..(n_uh2 + 1) {
                uh2_ordinates[(t - 1) as usize] = s_curves2(t as f64, self.x4) -
                    s_curves2(t as f64 - 1., self.x4);
            }

            for (p, e) in precip.iter().zip(potential_evap) {
                let net_evap;
                let mut routing_pattern;
                let reservoir_production;
                if p > e {
                    net_evap = 0.;
                    let mut scaled_net_precip = (p - e) / self.x1;
                    if scaled_net_precip > 13. {
                        scaled_net_precip = 13.;
                    }
                    let tanh_scaled_net_precip = scaled_net_precip.tanh();
                    reservoir_production = (self.x1 * (1. - (self.production_store / self.x1).powf(2.)) *
                                                tanh_scaled_net_precip) /
                        (1. + self.production_store / self.x1 * tanh_scaled_net_precip);

                    routing_pattern = p - e - reservoir_production;
                } else {
                    let mut scaled_net_evap = (e - p) / self.x1;
                    if scaled_net_evap > 13. {
                        scaled_net_evap = 13.;
                    }
                    let tanh_scaled_net_evap = scaled_net_evap.tanh();

                    let ps_div_x1 = (2. - self.production_store / self.x1) * tanh_scaled_net_evap;
                    net_evap = self.production_store * (ps_div_x1) /
                        (1. + (1. - self.production_store / self.x1) * tanh_scaled_net_evap);

                    reservoir_production = 0.;
                    routing_pattern = 0.;
                }

                self.production_store = self.production_store - net_evap + reservoir_production;

                let percolation = self.production_store /
                    (1. + (self.production_store / 2.25 / self.x1).powf(4.)).powf(0.25);
                routing_pattern = routing_pattern + (self.production_store - percolation);

                self.production_store = percolation;

                for i in 0..(self.uh1.len() - 1) {
                    self.uh1[i] = self.uh1[i + 1] + uh1_ordinates[i] * routing_pattern;
                }
                if let (Some(last_uh1), Some(last_ordinate)) = (self.uh1.last_mut(), uh1_ordinates.last()) {
                    *last_uh1 = *last_ordinate * routing_pattern;
                }

                for j in 0..(self.uh2.len() - 1) {
                    self.uh2[j] = self.uh2[j + 1] + uh2_ordinates[j] * routing_pattern
                }
                if let (Some(last_uh2), Some(last_ordinate)) = (self.uh2.last_mut(), uh2_ordinates.last()) {
                    *last_uh2 = *last_ordinate * routing_pattern;
                }

                let groundwater_exchange = self.x2 * (self.routing_store / self.x3).powf(3.5);
                self.routing_store = (0.0f64).max(self.routing_store + self.uh1[0] * 0.9 + groundwater_exchange);

                let r2 = self.routing_store / (1. + (self.routing_store / self.x3).powf(4.)).powf(0.25);
                let qr = self.routing_store - r2;
                self.routing_store = r2;
                let qd = (0.0f64).max(self.uh2[0] * 0.1 + groundwater_exchange);
                let q = qr + qd;

                self.qsim.push(q);
            }
        }

        pub fn init(
            &mut self,
            params: HashMap<&str, f64>,
            production_store: Option<f64>,
            routing_store: Option<f64>,
            unit_hydrographs: Option<HashMap<&str, Vec<f64>>>,
        ) {
            self.x1 = params["X1"];
            self.x2 = params["X2"];
            self.x3 = params["X3"];
            self.x4 = params["X4"];

            if production_store.is_some() {
                self.production_store = production_store.unwrap();
            }

            if routing_store.is_some() {
                self.routing_store = routing_store.unwrap();
            }

            let n_uh1 = self.x4.ceil() as usize;
            let n_uh2 = (2.0 * self.x4).ceil() as usize;

            // Default to all zeroes but set after if defined in unit_hydrographs HashMap
            self.uh1 = vec![0.;n_uh1];
            self.uh2 = vec![0.;n_uh2];

            if unit_hydrographs.is_some() {
                let unit_hydrographs = unit_hydrographs.unwrap();
                if unit_hydrographs.get("uh1").is_some() {
                    self.uh1.clone_from(&*unit_hydrographs.get("uh1").unwrap());
                }

                if unit_hydrographs.get("uh2").is_some() {
                    self.uh2.clone_from(&*unit_hydrographs.get("uh2").unwrap());
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use hydromodels;

    #[test]
    fn gr4j_simple_test() {
        let mut params = HashMap::new();

        params.insert("X1", 10.);
        params.insert("X2", 5.);
        params.insert("X3", 4.);
        params.insert("X4", 1.);

        let mut expected: Vec<f64> = Vec::new();
        expected.push(0.13030636843636356);
        expected.push(0.7348755690383941);
        expected.push(3.8482364547176102);
        expected.push(11.246676991863168);
        expected.push(13.90322269162079);

        let mut gr4j = hydromodels::GR4JModel {
            .. Default::default()
        };

        gr4j.init(params, None, None, None);
        gr4j.run(
            &[10., 2., 3., 4., 5.],
            &[0.5, 0.5, 0.5, 0.5, 0.5],
        );

        assert_eq!(gr4j.qsim, expected);
    }

    #[test]
    fn gr4j_scaled_test() {
        let mut params = HashMap::new();

        params.insert("X1", 10.);
        params.insert("X2", 5.);
        params.insert("X3", 4.);
        params.insert("X4", 1.);

        let mut expected: Vec<f64> = Vec::new();

        expected.push(0.13030636843636356);
        expected.push(0.5196067786804051);
        expected.push(0.7978026532390882);
        expected.push(128.2406932741725);
        expected.push(20.992238476264593);

        let mut gr4j = hydromodels::GR4JModel {
            .. Default::default()
        };

        gr4j.init(params, None, None, None);
        gr4j.run(
            &[10., 2., 3., 150., 5.],
            &[0.5, 14., 0.5, 10., 0.5],
        );

        assert_eq!(gr4j.qsim, expected);
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

        let mut expected: Vec<f64> = Vec::new();

        expected.push(5.1602235324393675);
        expected.push(10.091188499285725);
        expected.push(9.82974398339987);
        expected.push(136.95699908376093);
        expected.push(21.019904684254975);

        let mut gr4j = hydromodels::GR4JModel {
            .. Default::default()
        };

        gr4j.init(params, Some(10.), None, None);
        gr4j.run(
            &[10., 2., 3., 150., 5.],
            &[0.5, 14., 0.5, 10., 0.5],
        );

        assert_eq!(gr4j.qsim, expected);
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

        let mut expected: Vec<f64> = Vec::new();

        expected.push(7.031781387527497);
        expected.push(15.758863927257103);
        expected.push(10.450019503728232);
        expected.push(32.38841274927161);
        expected.push(115.58026737087125);

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
            .. Default::default()
        };

        gr4j.init(params, Some(10.), None, Some(unit_hydrographs));
        gr4j.run(
            &[10., 2., 3., 150., 5.],
            &[0.5, 14., 0.5, 10., 0.5],
        );

        assert_eq!(gr4j.qsim, expected);
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

        let mut expected: Vec<f64> = Vec::new();

        expected.push(1.326727594436605);
        expected.push(14.387525661905979);
        expected.push(10.409917456973545);
        expected.push(32.38841274927161);
        expected.push(115.58026737087125);

        let mut unit_hydrographs = HashMap::new();
        unit_hydrographs.insert("uh1", vec![111.13119599074196, 3.7349877368581]);

        let mut gr4j = hydromodels::GR4JModel {
            .. Default::default()
        };

        gr4j.init(params, Some(10.), None, Some(unit_hydrographs));
        gr4j.run(
            &[10., 2., 3., 150., 5.],
            &[0.5, 14., 0.5, 10., 0.5],
        );

        assert_eq!(gr4j.qsim, expected);
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
