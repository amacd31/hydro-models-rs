mod models;

pub use models::gr2m::{GR2MParams, GR2MModel};
pub use models::gr4j::{GR4JParams, GR4JModel};

#[cfg(test)]
mod tests {
    use crate::models::gr4j;
    use std::collections::HashMap;

    #[test]
    fn gr4j_create_test() {
        let params = crate::GR4JParams {
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

        let mut gr4j = crate::GR4JModel::create(params);

        let qsim = gr4j.run(&[10., 2., 3., 4., 5.], &[0.5, 0.5, 0.5, 0.5, 0.5]);

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

        let mut gr4j = crate::GR4JModel {
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

        let mut gr4j = crate::GR4JModel {
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

        let mut gr4j = crate::GR4JModel {
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

        let mut gr4j = crate::GR4JModel {
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

        let mut gr4j = crate::GR4JModel {
            ..Default::default()
        };

        gr4j.init(&params, Some(10.), None, Some(unit_hydrographs));
        let qsim = gr4j.run(&[10., 2., 3., 150., 5.], &[0.5, 14., 0.5, 10., 0.5]);

        assert_eq!(qsim, expected);
    }

    #[test]
    fn s_curves1_test() {
        assert_eq!(gr4j::s_curves1(-1.0, 2.0), 0.0);
        assert_eq!(gr4j::s_curves1(0.0, 2.0), 0.0);
        assert_eq!(gr4j::s_curves1(1.0, 2.0), 0.1767766952966369);
        assert_eq!(gr4j::s_curves1(2.0, 2.0), 1.0);
        assert_eq!(gr4j::s_curves1(2.1, 2.0), 1.0);
    }

    #[test]
    fn s_curves2_test() {
        assert_eq!(gr4j::s_curves2(-1.0, 2.0), 0.0);
        assert_eq!(gr4j::s_curves2(0.0, 2.0), 0.0);
        assert_eq!(gr4j::s_curves2(1.0, 2.0), 0.08838834764831845);
        assert_eq!(gr4j::s_curves2(2.0, 2.0), 0.5);
        assert_eq!(gr4j::s_curves2(2.1, 2.0), 0.5601759051904955);
        assert_eq!(gr4j::s_curves2(4.0, 2.0), 1.0);
        assert_eq!(gr4j::s_curves2(4.1, 2.0), 1.0);
    }

    #[test]
    fn gr2m_test() {
        let params = crate::GR2MParams { x1: 1., x2: 1. };

        let expected = vec![
            0.0009450053530675536,
            0.0327630181266177,
            0.2055971698811445,
            0.6632149029988136,
            0.991927371887451,
            1.0533436405957102,
            0.9029856114321761,
        ];

        let mut gr2m = crate::GR2MModel::create(params);

        let qsim = gr2m.run(&[1., 2., 3., 4., 3., 2., 1.], &[1., 1., 1., 1., 1., 1., 1.]);

        assert_eq!(gr2m.production_store, 0.1792526700466929);
        assert_eq!(gr2m.routing_store, 6.922989036389427);
        assert_eq!(qsim, expected);
    }
}
