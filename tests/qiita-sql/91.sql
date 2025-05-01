-- Adding 3 studies and their sample-metadata to test qp-knight-lab-processing pluing
INSERT INTO
    qiita.study (
        study_id,
        email,
        timeseries_type_id,
        metadata_complete,
        mixs_compliant,
        principal_investigator_id,
        reprocess,
        study_title,
        study_alias,
        study_description,
        study_abstract
    )
VALUES
    (
        6123,
        'demo@microbio.me',
        1,
        true,
        true,
        1,
        false,
        'Study 6123',
        '',
        '',
        ''
    ),
    (
        6124,
        'demo@microbio.me',
        1,
        true,
        true,
        1,
        false,
        'Study 6124',
        '',
        '',
        ''
    ),
    (
        11661,
        'demo@microbio.me',
        1,
        true,
        true,
        1,
        false,
        'Study 11661',
        '',
        '',
        ''
    ),
    (
        13059,
        'demo@microbio.me',
        1,
        true,
        true,
        1,
        false,
        'Study 13059',
        '',
        '',
        ''
    );

-- Add studies to Qiita portal
INSERT INTO qiita.study_portal (study_id, portal_type_id) VALUES (6123, 1), (6124, 1), (11661, 1), (13059, 1);

-- Creating the sample-metadata for tne new studies

CREATE TABLE qiita.sample_6123 (
    sample_id VARCHAR NOT NULL PRIMARY KEY,
    sample_values JSONB NOT NULL
);
INSERT INTO qiita.sample_6123 (sample_id, sample_values) VALUES  ('qiita_sample_column_names', '{}'::json);

CREATE TABLE qiita.sample_6124 (
    sample_id VARCHAR NOT NULL PRIMARY KEY,
    sample_values JSONB NOT NULL
);
INSERT INTO qiita.sample_6124 (sample_id, sample_values) VALUES  ('qiita_sample_column_names', '{}'::json);

CREATE TABLE qiita.sample_11661 (
    sample_id VARCHAR NOT NULL PRIMARY KEY,
    sample_values JSONB NOT NULL
);
INSERT INTO qiita.sample_11661 (sample_id, sample_values) VALUES  ('qiita_sample_column_names', '{}'::json);

CREATE TABLE qiita.sample_13059 (
    sample_id VARCHAR NOT NULL PRIMARY KEY,
    sample_values JSONB NOT NULL
);
INSERT INTO qiita.sample_13059 (sample_id, sample_values) VALUES  ('qiita_sample_column_names', '{}'::json);

-- Now adding the samples
INSERT INTO
    qiita.study_sample (study_id, sample_id)
VALUES
    -- 11661
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-143'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-144'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-145'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-146'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-147'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-148'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-149'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-150'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-151'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-152'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-153'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-154'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-155'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-156'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-157'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-158'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-159'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-160'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-161'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-162'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-163'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-164'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-165'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-166'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-167'),
    (11661, '11661.CDPH-SAL.Salmonella.Typhi.MDL-168'),
    (11661, '11661.P21.E.coli.ELI344'),
    (11661, '11661.P21.E.coli.ELI345'),
    (11661, '11661.P21.E.coli.ELI347'),
    (11661, '11661.P21.E.coli.ELI348'),
    (11661, '11661.P21.E.coli.ELI349'),
    (11661, '11661.P21.E.coli.ELI350'),
    (11661, '11661.P21.E.coli.ELI351'),
    (11661, '11661.P21.E.coli.ELI352'),
    (11661, '11661.P21.E.coli.ELI353'),
    (11661, '11661.P21.E.coli.ELI354'),
    (11661, '11661.P21.E.coli.ELI355'),
    (11661, '11661.P21.E.coli.ELI357'),
    (11661, '11661.P21.E.coli.ELI358'),
    (11661, '11661.P21.E.coli.ELI359'),
    (11661, '11661.P21.E.coli.ELI361'),
    (11661, '11661.P21.E.coli.ELI362'),
    (11661, '11661.P21.E.coli.ELI363'),
    (11661, '11661.P21.E.coli.ELI364'),
    (11661, '11661.P21.E.coli.ELI365'),
    (11661, '11661.P21.E.coli.ELI366'),
    (11661, '11661.P21.E.coli.ELI367'),
    (11661, '11661.P21.E.coli.ELI368'),
    (11661, '11661.P21.E.coli.ELI369'),
    (11661, '11661.stALE.E.coli.A1.F21.I1.R1'),
    (11661, '11661.stALE.E.coli.A2.F21.I1.R1'),
    (11661, '11661.stALE.E.coli.A3.F18.I1.R1'),
    (11661, '11661.stALE.E.coli.A3.F40.I1.R1'),
    (11661, '11661.stALE.E.coli.A4.F21.I1.R1'),
    (11661, '11661.stALE.E.coli.A4.F21.I1.R2'),
    (11661, '11661.stALE.E.coli.A4.F42.I1.R1'),
    (11661, '11661.stALE.E.coli.A5.F21.I1.R1'),
    (11661, '11661.stALE.E.coli.A5.F42.I1.R1'),
    (11661, '11661.stALE.E.coli.A6.F21.I1.R1'),
    (11661, '11661.stALE.E.coli.A6.F43.I1.R1'),
    (11661, '11661.stALE.E.coli.A7.F21.I1.R1'),
    (11661, '11661.stALE.E.coli.A7.F42.I1.R1'),
    (11661, '11661.stALE.E.coli.A8.F20.I1.R1'),
    (11661, '11661.stALE.E.coli.A8.F42.I1.R1'),
    (11661, '11661.stALE.E.coli.A9.F21.I1.R1'),
    (11661, '11661.stALE.E.coli.A9.F44.I1.R1'),
    (11661, '11661.stALE.E.coli.A10.F21.I1.R1'),
    (11661, '11661.stALE.E.coli.A10.F43.I1.R1'),
    (11661, '11661.stALE.E.coli.A10.F131.I1.R1'),
    (11661, '11661.stALE.E.coli.A11.F21.I1.R1'),
    (11661, '11661.stALE.E.coli.A11.F43.I1.R1'),
    (11661, '11661.stALE.E.coli.A11.F119.I1.R1'),
    (11661, '11661.stALE.E.coli.A12.F21.I1.R1'),
    (11661, '11661.stALE.E.coli.A12.F43.I1.R1'),
    (11661, '11661.stALE.E.coli.A12.F136.I1.R1'),
    (11661, '11661.stALE.E.coli.A13.F20.I1.R1'),
    (11661, '11661.stALE.E.coli.A13.F42.I1.R1'),
    (11661, '11661.stALE.E.coli.A13.F121.I1.R1'),
    (11661, '11661.stALE.E.coli.A14.F20.I1.R1'),
    (11661, '11661.stALE.E.coli.A14.F42.I1.R1'),
    (11661, '11661.stALE.E.coli.A14.F133.I1.R1'),
    (11661, '11661.stALE.E.coli.A15.F21.I1.R1'),
    (11661, '11661.stALE.E.coli.A15.F42.I1.R1'),
    (11661, '11661.stALE.E.coli.A15.F117.I1.R1'),
    (11661, '11661.stALE.E.coli.A16.F20.I1.R1'),
    (11661, '11661.stALE.E.coli.A16.F42.I1.R1'),
    (11661, '11661.stALE.E.coli.A16.F134.I1.R1'),
    (11661, '11661.stALE.E.coli.A17.F21.I1.R1'),
    (11661, '11661.stALE.E.coli.A17.F118.I1.R1'),
    (11661, '11661.stALE.E.coli.A18.F18.I1.R1'),
    (11661, '11661.stALE.E.coli.A18.F39.I1.R1'),
    (11661, '11661.stALE.E.coli.A18.F130.I1.R1'),
    (11661, '11661.BLANK.40.12G'),
    (11661, '11661.BLANK.40.12H'),
    (11661, '11661.Pputida.JBEI.HGL.Pputida.107.BP6'),
    (11661, '11661.Pputida.JBEI.HGL.Pputida.108.BP7'),
    (11661, '11661.Pputida.JBEI.HGL.Pputida.109.BP8'),
    (11661, '11661.Pputida.JBEI.HGL.Pputida.110.M2'),
    (11661, '11661.Pputida.JBEI.HGL.Pputida.111.M5'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.112'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.113'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.114'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.115'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.116'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.117'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.118'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.119'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.120'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.121'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.122'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.123'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.124'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.125'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.126'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.127'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.128'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.129'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.130'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.131'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.132'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.133'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.134'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.135'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.136'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.137'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.138'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.139'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.140'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.141'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.142'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.143'),
    (11661, '11661.Pputida.TALE.HGL.Pputida.144'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.145'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.146'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.147'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.148'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.149'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.150'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.151'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.152'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.153'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.154'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.155'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.156'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.157'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.158'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.159'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.160'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.161'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.162'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.163'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.164'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.165'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.166'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.167'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.168'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.169'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.170'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.171'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.172'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.173'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.174'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.175'),
    (11661, '11661.Pputida.PALE.HGL.Pputida.176'),
    (11661, '11661.JM-Metabolic.GN0.2005'),
    (11661, '11661.JM-Metabolic.GN0.2007'),
    (11661, '11661.JM-Metabolic.GN0.2009'),
    (11661, '11661.JM-Metabolic.GN0.2094'),
    (11661, '11661.JM-Metabolic.GN0.2099'),
    (11661, '11661.JM-Metabolic.GN0.2148'),
    (11661, '11661.JM-Metabolic.GN0.2165'),
    (11661, '11661.JM-Metabolic.GN0.2169'),
    (11661, '11661.JM-Metabolic.GN0.2172'),
    (11661, '11661.JM-Metabolic.GN0.2175'),
    (11661, '11661.JM-Metabolic.GN0.2183'),
    (11661, '11661.JM-Metabolic.GN0.2215'),
    (11661, '11661.JM-Metabolic.GN0.2254'),
    (11661, '11661.JM-Metabolic.GN0.2277'),
    (11661, '11661.JM-Metabolic.GN0.2290'),
    (11661, '11661.JM-Metabolic.GN0.2337'),
    (11661, '11661.JM-Metabolic.GN0.2317'),
    (11661, '11661.JM-Metabolic.GN0.2354'),
    (11661, '11661.JM-Metabolic.GN0.2375'),
    (11661, '11661.JM-Metabolic.GN0.2380'),
    (11661, '11661.JM-Metabolic.GN0.2393'),
    (11661, '11661.JM-Metabolic.GN0.2404'),
    (11661, '11661.BLANK.41.12H'),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.BOP27.4.14'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.BOP27.4.23'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.BOP27.4.48'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.BOP27.6.21'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.BOP27.6.35'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.BOP27.10.13'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.BOP27.10.28'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.BOP27.10.51'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.Lib4.18.19'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.Lib4.18.59'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.Lib4.18.35'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.Lib4.20.16'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.Lib4.20.43'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.Lib4.20.71'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.Lib4.22.16'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.Lib4.22.28'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.Lib4.22.52'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.Lib4.24.9'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.Lib4.24.24'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.Lib4.24.52'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.Lib4.26.6'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.Lib4.26.27'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.Lib4.26.69'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.Lib4.28.13'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.Lib4.28.28'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.Lib4.28.53'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.Lib4.30.7'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.Lib4.30.22'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.Lib4.30.60'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.Lib4.32.6'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.Lib4.32.20'
    ),
    (
        11661,
        '11661.Deoxyribose.PALE.ALE.MG1655.Lib4.32.56'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.1.24'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.1.57'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.1.69'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.3.23'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.3.50'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.3.61'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.5.22'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.5.36'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.5.46'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.7.23'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.7.41'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.7.51'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.17.25'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.17.58'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.17.64'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.19.25'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.19.55'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.19.63'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.21.23'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.21.46'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.21.51'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.29.25'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.29.49'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.29.57'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.31.24'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.31.42'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.31.62'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.33.21'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.33.41'
    ),
    (
        11661,
        '11661.AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.33.50'
    ),
    (11661, '11661.JM-Metabolic.GN02514'),
    (11661, '11661.JM-Metabolic.GN02529'),
    (11661, '11661.JM-Metabolic.GN02531'),
    (11661, '11661.JM-Metabolic.GN02567'),
    (11661, '11661.JM-Metabolic.GN02590'),
    (11661, '11661.JM-Metabolic.GN02657'),
    (11661, '11661.JM-Metabolic.GN02748'),
    (11661, '11661.JM-Metabolic.GN02766'),
    (11661, '11661.JM-Metabolic.GN02769'),
    (11661, '11661.JM-Metabolic.GN02787'),
    (11661, '11661.JM-Metabolic.GN03132'),
    (11661, '11661.JM-Metabolic.GN03218'),
    (11661, '11661.JM-Metabolic.GN03252'),
    (11661, '11661.JM-Metabolic.GN03409'),
    (11661, '11661.JM-Metabolic.GN04014'),
    (11661, '11661.JM-Metabolic.GN04094'),
    (11661, '11661.JM-Metabolic.GN04255'),
    (11661, '11661.JM-Metabolic.GN04306'),
    (11661, '11661.JM-Metabolic.GN04428'),
    (11661, '11661.JM-Metabolic.GN04488'),
    (11661, '11661.JM-Metabolic.GN04540'),
    (11661, '11661.JM-Metabolic.GN04563'),
    (11661, '11661.JM-Metabolic.GN04612'),
    (11661, '11661.JM-Metabolic.GN04665'),
    (11661, '11661.JM-Metabolic.GN04682'),
    (11661, '11661.JM-Metabolic.GN05002'),
    (11661, '11661.JM-Metabolic.GN05109'),
    (11661, '11661.JM-Metabolic.GN05128'),
    (11661, '11661.JM-Metabolic.GN05367'),
    (11661, '11661.JM-Metabolic.GN05377'),
    (11661, '11661.BLANK.42.12G'),
    (11661, '11661.BLANK.42.12H'),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0326'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0327'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0328'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0329'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0330'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0352'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0353'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0354'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0355'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0356'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0357'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0364'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0366'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0367'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0368'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0369'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0370'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0371'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0372'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0373'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0374'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0375'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0376'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0377'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0378'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0380'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0381'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0382'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0383'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0384'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0385'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0386'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0387'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0388'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0389'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0390'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0391'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0392'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0393'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0394'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0395'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0396'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0397'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0398'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0399'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0400'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0401'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0402'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0403'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0404'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0405'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0406'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0407'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0408'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0409'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0417'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0418'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0419'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0420'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0421'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0473'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0474'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0483'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0484'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0485'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0486'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0516'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0517'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0518'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0519'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0520'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0521'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0522'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0523'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0524'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-B0525'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-R08624'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-R08704'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-R10727'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-R11044'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-R11078'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-R11101'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-R11102'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-R11103'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-R11135'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-R11153'
    ),
    (
        11661,
        '11661.JM-MEC.Staphylococcus.aureusstrain.BERTI-R11154'
    ),
    (11661, '11661.JM-Metabolic.GN02424'),
    (11661, '11661.JM-Metabolic.GN02446'),
    (11661, '11661.JM-Metabolic.GN02449'),
    (11661, '11661.JM-Metabolic.GN02487'),
    (11661, '11661.JM-Metabolic.GN02501'),
    (11661, '11661.BLANK.43.12G'),
    (11661, '11661.BLANK.43.12H'),
    (11661, '11661.RMA.KHP.rpoS.Mage.Q97D'),
    (11661, '11661.RMA.KHP.rpoS.Mage.Q97L'),
    (11661, '11661.RMA.KHP.rpoS.Mage.Q97N'),
    (11661, '11661.RMA.KHP.rpoS.Mage.Q97E'),
    (11661, '11661.JBI.KHP.HGL.021'),
    (11661, '11661.JBI.KHP.HGL.022'),
    (11661, '11661.JBI.KHP.HGL.023'),
    (11661, '11661.JBI.KHP.HGL.024'),
    (11661, '11661.JBI.KHP.HGL.025'),
    (11661, '11661.JBI.KHP.HGL.026'),
    (11661, '11661.JBI.KHP.HGL.027'),
    (11661, '11661.JBI.KHP.HGL.028.Amitesh.soxR'),
    (11661, '11661.JBI.KHP.HGL.029.Amitesh.oxyR'),
    (11661, '11661.JBI.KHP.HGL.030.Amitesh.soxR.oxyR'),
    (11661, '11661.JBI.KHP.HGL.031.Amitesh.rpoS'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_Lib4_24_52'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_128'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-150'),
    (11661, '11661.JM-Metabolic__GN02531'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_31_62'),
    (11661, '11661.stALE_E_coli_A14_F42_I1R1'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_BOP27_4_14'),
    (11661, '11661.JM-Metabolic__GN0_2254'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_1_69'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_Lib4_24_24'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_147'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-162'),
    (11661, '11661.P21_E_coli_ELI364'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_19_63'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_155'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_124'),
    (11661, '11661.JM-Metabolic__GN05109'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_173'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_Lib4_26_69'),
    (11661, '11661.P21_E_coli_ELI347'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-156'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_33_41'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_Lib4_32_20'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_Lib4_22_16'),
    (11661, '11661.JM-Metabolic__GN04540'),
    (11661, '11661.stALE_E_coli_A5_F21_I1R1'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_33_21'),
    (11661, '11661.P21_E_coli_ELI368'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_151'),
    (11661, '11661.P21_E_coli_ELI349'),
    (11661, '11661.stALE_E_coli_A17_F118_I1R1'),
    (11661, '11661.stALE_E_coli_A12_F43_I1R1'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_BOP27_10_13'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_160'),
    (11661, '11661.JM-Metabolic__GN0_2337'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_7_41'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-145'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-159'),
    (11661, '11661.P21_E_coli_ELI355'),
    (11661, '11661.P21_E_coli_ELI353'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_166'),
    (11661, '11661.P21_E_coli_ELI369'),
    (11661, '11661.JM-Metabolic__GN0_2007'),
    (11661, '11661.Pputida_JBEI__HGL_Pputida_111_M5'),
    (11661, '11661.JM-Metabolic__GN02657'),
    (11661, '11661.JM-Metabolic__GN0_2290'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_136'),
    (11661, '11661.Pputida_JBEI__HGL_Pputida_108_BP7'),
    (11661, '11661.stALE_E_coli_A15_F117_I1R1'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_7_23'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_131'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_Lib4_18_35'),
    (11661, '11661.P21_E_coli_ELI366'),
    (11661, '11661.P21_E_coli_ELI363'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_148'),
    (11661, '11661.JM-Metabolic__GN04665'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-147'),
    (11661, '11661.JM-Metabolic__GN0_2172'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_Lib4_22_28'),
    (11661, '11661.P21_E_coli_ELI344'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_122'),
    (11661, '11661.stALE_E_coli_A5_F42_I1R1'),
    (11661, '11661.P21_E_coli_ELI357'),
    (11661, '11661.stALE_E_coli_A13_F121_I1R1'),
    (11661, '11661.JM-Metabolic__GN0_2215'),
    (11661, '11661.JM-Metabolic__GN05367'),
    (11661, '11661.JM-Metabolic__GN0_2354'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_156'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-163'),
    (11661, '11661.JM-Metabolic__GN05002'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_170'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_7_51'),
    (11661, '11661.JM-Metabolic__GN0_2165'),
    (11661, '11661.Pputida_JBEI__HGL_Pputida_110_M2'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_21_46'),
    (11661, '11661.JM-Metabolic__GN0_2005'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_BOP27_6_35'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_21_23'),
    (11661, '11661.P21_E_coli_ELI358'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_19_25'),
    (11661, '11661.P21_E_coli_ELI365'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_Lib4_18_19'),
    (11661, '11661.JM-Metabolic__GN03409'),
    (11661, '11661.JM-Metabolic__GN03252'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-148'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_Lib4_22_52'),
    (11661, '11661.JM-Metabolic__GN0_2380'),
    (11661, '11661.stALE_E_coli_A4_F42_I1R1'),
    (11661, '11661.P21_E_coli_ELI359'),
    (11661, '11661.JM-Metabolic__GN0_2009'),
    (11661, '11661.stALE_E_coli_A6_F43_I1R1'),
    (11661, '11661.JM-Metabolic__GN03132'),
    (11661, '11661.stALE_E_coli_A10_F43_I1R1'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_31_42'),
    (11661, '11661.P21_E_coli_ELI362'),
    (11661, '11661.P21_E_coli_ELI350'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_126'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_157'),
    (11661, '11661.stALE_E_coli_A15_F21_I1R1'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_133'),
    (11661, '11661.P21_E_coli_ELI361'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_164'),
    (11661, '11661.stALE_E_coli_A4_F21_I1R1'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_BOP27_4_23'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-160'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_Lib4_20_16'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-167'),
    (11661, '11661.stALE_E_coli_A12_F21_I1R1'),
    (11661, '11661.stALE_E_coli_A12_F136_I1R1'),
    (11661, '11661.JM-Metabolic__GN04682'),
    (11661, '11661.JM-Metabolic__GN04306'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_154'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-168'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_112'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_1_57'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_175'),
    (11661, '11661.JM-Metabolic__GN0_2317'),
    (11661, '11661.stALE_E_coli_A3_F40_I1R1'),
    (11661, '11661.stALE_E_coli_A16_F134_I1R1'),
    (11661, '11661.JM-Metabolic__GN05377'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-164'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_119'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_159'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_118'),
    (11661, '11661.stALE_E_coli_A7_F21_I1R1'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_169'),
    (11661, '11661.stALE_E_coli_A18_F39_I1R1'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-166'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-153'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_29_57'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-158'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_145'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_Lib4_20_43'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_19_55'),
    (11661, '11661.stALE_E_coli_A11_F21_I1R1'),
    (11661, '11661.stALE_E_coli_A3_F18_I1R1'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_Lib4_24_9'),
    (11661, '11661.P21_E_coli_ELI348'),
    (11661, '11661.JM-Metabolic__GN0_2277'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_Lib4_28_53'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_21_51'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-165'),
    (11661, '11661.stALE_E_coli_A10_F21_I1R1'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_138'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-151'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_161'),
    (11661, '11661.stALE_E_coli_A8_F42_I1R1'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_140'),
    (11661, '11661.stALE_E_coli_A14_F133_I1R1'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_132'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-144'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_121'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_129'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_120'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_168'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_165'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_BOP27_10_51'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_3_23'),
    (11661, '11661.stALE_E_coli_A16_F42_I1R1'),
    (11661, '11661.Pputida_JBEI__HGL_Pputida_109_BP8'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_Lib4_28_28'),
    (11661, '11661.P21_E_coli_ELI367'),
    (11661, '11661.stALE_E_coli_A14_F20_I1R1'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_142'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_3_61'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_130'),
    (11661, '11661.stALE_E_coli_A1_F21_I1R1'),
    (11661, '11661.JM-Metabolic__GN04563'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_114'),
    (11661, '11661.JM-Metabolic__GN0_2393'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-146'),
    (11661, '11661.stALE_E_coli_A16_F20_I1R1'),
    (11661, '11661.stALE_E_coli_A18_F130_I1R1'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_Lib4_30_7'),
    (11661, '11661.JM-Metabolic__GN03218'),
    (11661, '11661.JM-Metabolic__GN02567'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_123'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_Lib4_26_6'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_150'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-157'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_33_50'),
    (11661, '11661.JM-Metabolic__GN02766'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_64'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_5_22'),
    (11661, '11661.P21_E_coli_ELI345'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_135'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_31_24'),
    (11661, '11661.stALE_E_coli_A9_F44_I1R1'),
    (11661, '11661.JM-Metabolic__GN02514'),
    (11661, '11661.JM-Metabolic__GN04488'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_58'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_115'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_174'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_139'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_143'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_134'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_Lib4_26_27'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-154'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_158'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-155'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_25'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_152'),
    (11661, '11661.stALE_E_coli_A4_F21_I1R2'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_171'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_153'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_29_49'),
    (11661, '11661.JM-Metabolic__GN02590'),
    (11661, '11661.JM-Metabolic__GN0_2169'),
    (11661, '11661.JM-Metabolic__GN05128'),
    (11661, '11661.stALE_E_coli_A6_F21_I1R1'),
    (11661, '11661.stALE_E_coli_A9_F21_I1R1'),
    (11661, '11661.JM-Metabolic__GN02769'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_BOP27_6_21'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_146'),
    (11661, '11661.P21_E_coli_ELI352'),
    (11661, '11661.stALE_E_coli_A18_F18_I1R1'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_Lib4_30_22'),
    (11661, '11661.JM-Metabolic__GN0_2148'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_137'),
    (11661, '11661.stALE_E_coli_A15_F42_I1R1'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-161'),
    (11661, '11661.JM-Metabolic__GN04428'),
    (11661, '11661.JM-Metabolic__GN0_2183'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_163'),
    (11661, '11661.JM-Metabolic__GN04255'),
    (11661, '11661.stALE_E_coli_A13_F20_I1R1'),
    (11661, '11661.JM-Metabolic__GN0_2094'),
    (11661, '11661.JM-Metabolic__GN02748'),
    (11661, '11661.P21_E_coli_ELI351'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_Lib4_30_60'),
    (11661, '11661.stALE_E_coli_A8_F20_I1R1'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_Lib4_32_6'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_176'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_Lib4_20_71'),
    (11661, '11661.P21_E_coli_ELI354'),
    (11661, '11661.JM-Metabolic__GN0_2175'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_Lib4_32_56'),
    (11661, '11661.JM-Metabolic__GN02529'),
    (11661, '11661.JM-Metabolic__GN04014'),
    (11661, '11661.Pputida_JBEI__HGL_Pputida_107_BP6'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_172'),
    (11661, '11661.stALE_E_coli_A17_F21_I1R1'),
    (11661, '11661.stALE_E_coli_A2_F21_I1R1'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_144'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_113'),
    (11661, '11661.JM-Metabolic__GN04094'),
    (11661, '11661.stALE_E_coli_A11_F119_I1R1'),
    (11661, '11661.stALE_E_coli_A10_F131_I1R1'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_5_46'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_Lib4_28_13'),
    (11661, '11661.stALE_E_coli_A13_F42_I1R1'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_167'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_29_25'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-143'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_162'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_117'),
    (11661, '11661.JM-Metabolic__GN0_2404'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_Lib4_18_59'),
    (11661, '11661.JM-Metabolic__GN0_2099'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-149'),
    (11661, '11661.JM-Metabolic__GN0_2375'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_5_36'),
    (11661, '11661.JM-Metabolic__GN04612'),
    (11661, '11661.CDPH-SAL_Salmonella_Typhi_MDL-152'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_125'),
    (11661, '11661.stALE_E_coli_A11_F43_I1R1'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_141'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_127'),
    (11661, '11661.Pputida_PALE__HGL_Pputida_149'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_BOP27_10_28'),
    (11661, '11661.Deoxyribose_PALE_ALE__MG1655_BOP27_4_48'),
    (11661, '11661.Pputida_TALE__HGL_Pputida_116'),
    (11661, '11661.stALE_E_coli_A7_F42_I1R1'),
    (11661, '11661.JM-Metabolic__GN02787'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_3_50'),
    (11661, '11661.AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_1_24'),

    -- 6123
    (6123, '6123.3A'),
    (6123, '6123.4A'),
    (6123, '6123.5B'),
    (6123, '6123.6A'),
    (6123, '6123.BLANK.41.12G'),
    (6123, '6123.7A'),
    (6123, '6123.8A'),
    (6123, '6123.ISB'),
    (6123, '6123.GFR'),
    -- 13059
    (13059, '13059.EP808104A01'),
    (13059, '13059.AP173299B04'),
    (13059, '13059.EP587478B04'),
    (13059, '13059.ep256643b01'),
    (13059, '13059.EP455759B04'),
    (13059, '13059.EP606652B04'),
    (13059, '13059.AP046327B02'),
    (13059, '13059.EP316863B03'),
    (13059, '13059.EP061002B01'),
    (13059, '13059.EP786631A04'),
    (13059, '13059.SP205732A02'),
    (13059, '13059.SP232079A01'),
    (13059, '13059.EP554513B02'),
    (13059, '13059.LP128543A01'),
    (13059, '13059.SP246941A01'),
    (13059, '13059.EP182060B03'),
    (13059, '13059.AP967057A04'),
    (13059, '13059.EP320438B01'),
    (13059, '13059.SP506933A04'),
    (13059, '13059.SP515763A04'),
    (13059, '13059.SP612495A04'),
    (13059, '13059.EP159695B01'),
    (13059, '13059.SP317293A02'),
    (13059, '13059.SP404412A02'),
    (13059, '13059.SP232311A04'),
    (13059, '13059.AP911328B01'),
    (13059, '13059.EP636802A.1'),
    (13059, '13059.EP339057B02'),
    (13059, '13059.SP415023A02'),
    (13059, '13059.EP573313B01'),
    (13059, '13059.EP379938B01'),
    (13059, '13059.SP511289A02'),
    (13059, '13059.AP298002B02'),
    (13059, '13059.AP780167B02'),
    (13059, '13059.EP790021A04'),
    (13059, '13059.AP029018B01'),
    (13059, '13059.AP771472A04'),
    (13059, '13059.EP128904B02'),
    (13059, '13059.EP446610B02'),
    (13059, '13059.SP503615A02'),
    (13059, '13059.EP032410B02'),
    (13059, '13059.EP808110A04'),
    (13059, '13059.EP533426B03'),
    (13059, '13059.AP668631B04'),
    (13059, '13059.EP184255B04'),
    (13059, '13059.EP768164A02'),
    (13059, '13059.LP169879A01'),
    (13059, '13059.SP584551A08'),
    (13059, '13059.AP324642B04'),
    (13059, '13059.SP404403A02'),
    (13059, '13059.SP235189A01'),
    (13059, '13059.EP927462A02'),
    (13059, '13059.EP207041B01'),
    (13059, '13059.EP001624B01'),
    (13059, '13059.SP284096A02'),
    (13059, '13059.EP805337A01'),
    (13059, '13059.EP657260A01'),
    (13059, '13059.EP721390A04'),
    (13059, '13059.EP587476B04'),
    (13059, '13059.AP173301B04'),
    (13059, '13059.EP729434A01'),
    (13059, '13059.AP549678B01'),
    (13059, '13059.EP649737A03'),
    (13059, '13059.EP649653A04'),
    (13059, '13059.EP410041B01'),
    (13059, '13059.EP533429B04'),
    (13059, '13059.EP337325B04'),
    (13059, '13059.EP216516B04'),
    (13059, '13059.EP573310B01'),
    (13059, '13059.AP046324B02'),
    (13059, '13059.EP542578B.4'),
    (13059, '13059.SP365864A04'),
    (13059, '13059.LP128540A01'),
    (13059, '13059.EP996831B04'),
    (13059, '13059.AP062219B03'),
    (13059, '13059.EP447927B04'),
    (13059, '13059.EP554515B04'),
    (13059, '13059.AP309872B03'),
    (13059, '13059.AP787247B04'),
    (13059, '13059.EP533388B01'),
    (13059, '13059.EP190307B01'),
    (13059, '13059.LP154986A01'),
    (13059, '13059.EP451428B04'),
    (13059, '13059.LP191039A01'),
    (13059, '13059.EP282107B01'),
    (13059, '13059.EP905975A04'),
    (13059, '13059.EP981129A02'),
    (13059, '13059.LP166715A01'),
    (13059, '13059.EP617441B01'),
    (13059, '13059.EP927458A04'),
    (13059, '13059.EP090129B04'),
    (13059, '13059.EP587477B04'),
    (13059, '13059.SP230382A04'),
    (13059, '13059.SP229387A04'),
    (13059, '13059.EP657386A01'),
    (13059, '13059.EP649418A02'),
    (13059, '13059.EP260544B04'),
    (13059, '13059.SP640978A02'),
    (13059, '13059.SP230381A01'),
    (13059, '13059.AP744361A02'),
    (13059, '13059.SP584547A02'),
    (13059, '13059.AP696363B02'),
    (13059, '13059.EP987683A01'),
    (13059, '13059.EP163771B01'),
    (13059, '13059.EP333541B04'),
    (13059, '13059.SP573823A04'),
    (13059, '13059.EP606656B03'),
    (13059, '13059.SP232309A01'),
    (13059, '13059.EP718688A01'),
    (13059, '13059.EP940013A01'),
    (13059, '13059.EP899038A04'),
    (13059, '13059.SP415025A01'),
    (13059, '13059.EP667743A04'),
    (13059, '13059.SP577399A02'),
    (13059, '13059.EP617440B01'),
    (13059, '13059.EP230245B01'),
    (13059, '13059.EP784608A01'),
    (13059, '13059.EP970005A01'),
    (13059, '13059.EP484973B04'),
    (13059, '13059.EP685640B01'),
    (13059, '13059.EP032412B02'),
    (13059, '13059.AP470859B01'),
    (13059, '13059.SP612496A01'),
    (13059, '13059.EP447929B04'),
    (13059, '13059.EP447940B04'),
    (13059, '13059.C20'),
    (13059, '13059.EP256645B01'),
    (13059, '13059.EP479270B03'),
    (13059, '13059.SP754514A04'),
    (13059, '13059.LP128539A01'),
    (13059, '13059.EP980752B04'),
    (13059, '13059.EP929277A02'),
    (13059, '13059.EP749735A07'),
    (13059, '13059.EP431575B01'),
    (13059, '13059.LP128476A01'),
    (13059, '13059.EP724905B01'),
    (13059, '13059.EP291980B04'),
    (13059, '13059.EP504030B04'),
    (13059, '13059.EP790023A01'),
    (13059, '13059.C14'),
    (13059, '13059.EP808111A03'),
    (13059, '13059.EP400447B04'),
    (13059, '13059.AP891020A04'),
    (13059, '13059.SP631994A04'),
    (13059, '13059.EP970001A01'),
    (13059, '13059.22.001.710.503.791.00'),
    (13059, '13059.EP675042B01'),
    (13059, '13059.EP339061B02'),
    (13059, '13059.AP032413B04'),
    (13059, '13059.AP959450A03'),
    (13059, '13059.EP872341A01'),
    (13059, '13059.SP232310A04'),
    (13059, '13059.AP568785B04'),
    (13059, '13059.EP372981B04'),
    (13059, '13059.C18'),
    (13059, '13059.EP606662B04'),
    (13059, '13059.EP023801B04'),
    (13059, '13059.SP681591A04'),
    (13059, '13059.EP282276B04'),
    (13059, '13059.EP886422A01'),
    (13059, '13059.C9'),
    (13059, '13059.AP953594A02'),
    (13059, '13059.EP683835A01'),
    (13059, '13059.EP112567B02'),
    (13059, '13059.AP616837B04'),
    (13059, '13059.EP393712B02'),
    (13059, '13059.EP410042B01'),
    (13059, '13059.EP944059A02'),
    (13059, '13059.EP890157A02'),
    (13059, '13059.EP727972A04'),
    (13059, '13059.EP808106A01'),
    (13059, '13059.EP768748A04'),
    (13059, '13059.SP573824A04'),
    (13059, '13059.EP291979B04'),
    (13059, '13059.EP244366B01'),
    (13059, '13059.EP212214B01'),
    (13059, '13059.SP491898A02'),
    (13059, '13059.SP399724A04'),
    (13059, '13059.AP032412B04'),
    (13059, '13059.SP280481A02'),
    (13059, '13059.SP404405A02'),
    (13059, '13059.EP447926B04'),
    (13059, '13059.AP549681B02'),
    (13059, '13059.EP718687A04'),
    (13059, '13059.EP339059B02'),
    (13059, '13059.EP479266B04'),
    (13059, '13059.EP054632B01'),
    (13059, '13059.SP573859A04'),
    (13059, '13059.SP317297A02'),
    (13059, '13059.EP554518B04'),
    (13059, '13059.EP649623A01'),
    (13059, '13059.SP232270A02'),
    (13059, '13059.EP657385A04'),
    (13059, '13059.SP205754A01'),
    (13059, '13059.EP182065B04'),
    (13059, '13059.SP415021A02'),
    (13059, '13059.SP231629A02'),
    (13059, '13059.EP808105A01'),
    (13059, '13059.SP683466A02'),
    (13059, '13059.EP915769A04'),
    (13059, '13059.EP202095B04'),
    (13059, '13059.EP617442B01'),
    (13059, '13059.EP012991B03'),
    (13059, '13059.EP339053B02'),
    (13059, '13059.EP393714B01'),
    (13059, '13059.EP238034B01'),
    (13059, '13059.AP581451B02'),
    (13059, '13059.EP606663B04'),
    (13059, '13059.AP732307B04'),
    (13059, '13059.SP491907A02'),
    (13059, '13059.LP196272A01'),
    (13059, '13059.EP448041B04'),
    (13059, '13059.C3'),
    (13059, '13059.EP587475B04'),
    (13059, '13059.AP687591B04'),
    (13059, '13059.EP479894B04'),
    (13059, '13059.SP331134A04'),
    (13059, '13059.EP584756B04'),
    (13059, '13059.EP790019A01'),
    (13059, '13059.SP515443A04'),
    (13059, '13059.SP478193A02'),
    (13059, '13059.EP431562B04'),
    (13059, '13059.SP491897A02'),
    (13059, '13059.SP464350A04'),
    (13059, '13059.EP843906A04'),
    (13059, '13059.SP231628A02'),
    (13059, '13059.EP483291B04'),
    (13059, '13059.SP704319A04'),
    (13059, '13059.LP127767A01'),
    (13059, '13059.AP103463B01'),
    (13059, '13059.SP531696A04'),
    (13059, '13059.SP491900A02'),
    (13059, '13059.SP284095A03'),
    (13059, '13059.EP393718B01'),
    (13059, '13059.EP738469A01'),
    (13059, '13059.SP453872A01'),
    (13059, '13059.SP388683A02'),
    (13059, '13059.EP890158A02'),
    (13059, '13059.EP282108B01'),
    (13059, '13059.EP927461A04'),
    (13059, '13059.SP404409A02'),
    (13059, '13059.EP529635B02'),
    (13059, '13059.EP533389B03'),
    (13059, '13059.EP542577B04'),
    (13059, '13059.AP795068B04'),
    (13059, '13059.SP231630A02'),
    (13059, '13059.EP207042B04'),
    (13059, '13059.EP422407B01'),
    (13059, '13059.SP232077A04'),
    (13059, '13059.EP554506B04'),
    (13059, '13059.EP675044A01'),
    (13059, '13059.SP408629A01'),
    (13059, '13059.SP232114A04'),
    (13059, '13059.SP353893A02'),
    (13059, '13059.EP393717B01'),
    (13059, '13059.SP235186A04'),
    (13059, '13059.EP759450A04'),
    (13059, '13059.SP335002A04'),
    (13059, '13059.EP273332B04'),
    (13059, '13059.LP128541A01'),
    (13059, '13059.SP416130A04'),
    (13059, '13059.EP554501B04'),
    (13059, '13059.22.001.801.552.503.00'),
    (13059, '13059.EP446604B03'),
    (13059, '13059.EP121011B.1'),
    (13059, '13059.AP173305B04'),
    (13059, '13059.EP446602B.1'),
    (13059, '13059.SP231631A02'),
    (13059, '13059.C5'),
    (13059, '13059.EP656055A04'),
    (13059, '13059.EP702221B04'),
    (13059, '13059.EP128910B01'),
    (13059, '13059.SP490298A02'),
    (13059, '13059.EP868682A01'),
    (13059, '13059.LP154981A01'),
    (13059, '13059.EP455757B04'),
    (13059, '13059.SP471496A04'),
    (13059, '13059.AP568787B02'),
    (13059, '13059.SP573860A01'),
    (13059, '13059.SP257517A04'),
    (13059, '13059.EP073160B01'),
    (13059, '13059.EP772145A02'),
    (13059, '13059.EP073216B01'),
    (13059, '13059.EP617443B01'),
    (13059, '13059.SP573843A04'),
    (13059, '13059.EP385384B01'),
    (13059, '13059.AP006367B02'),
    (13059, '13059.EP790020A02'),
    (13059, '13059.EP487995B04'),
    (13059, '13059.SP415030A01'),
    (13059, '13059.EP393715B01'),
    (13059, '13059.EP410046B01'),
    (13059, '13059.AP745799A04'),
    (13059, '13059.AP668628B04'),
    (13059, '13059.AP223470B01'),
    (13059, '13059.EP921593A04'),
    (13059, '13059.SP645141A03'),
    (13059, '13059.SP641029A02'),
    (13059, '13059.LP127890A01'),
    (13059, '13059.EP337425B01'),
    (13059, '13059.EP182243B02'),
    (13059, '13059.EP073209B02'),
    (13059, '13059.EP207036B01'),
    (13059, '13059.EP729433A02'),
    (13059, '13059.EP479794B02'),
    (13059, '13059.EP447975B02'),
    (13059, '13059.AP065292B01'),
    (13059, '13059.SP257519A04'),
    (13059, '13059.EP846485A01'),
    (13059, '13059.EP808109A01'),
    (13059, '13059.LP128479A01'),
    (13059, '13059.EP921594A04'),
    (13059, '13059.EP738468A01'),
    (13059, '13059.EP675075A04'),
    (13059, '13059.EP447928B04'),
    (13059, '13059.C6'),
    (13059, '13059.EP023808B02'),
    (13059, '13059.EP244360B01'),
    (13059, '13059.EP573296B01'),
    (13059, '13059.EP256644B01'),
    (13059, '13059.EP182346B04'),
    (13059, '13059.EP772143A02'),
    (13059, '13059.EP385379B01'),
    (13059, '13059.SP230380A02'),
    (13059, '13059.EP876243A04'),
    (13059, '13059.EP260543B04'),
    (13059, '13059.EP455763B04'),
    (13059, '13059.EP305735B04'),
    (13059, '13059.EP808112A04'),
    (13059, '13059.SP511294A04'),
    (13059, '13059.SP247340A04'),
    (13059, '13059.lp127896a01'),
    (13059, '13059.EP882752A01'),
    (13059, '13059.EP431570B01'),
    (13059, '13059.AP905750A02'),
    (13059, '13059.SP561451A04'),
    (13059, '13059.AP481403B02'),
    (13059, '13059.EP202452B01'),
    (13059, '13059.EP001625B01'),
    (13059, '13059.LP128538A01'),
    (13059, '13059.EP927459A04'),
    (13059, '13059.EP159692B04'),
    (13059, '13059.SP410793A01'),
    (13059, '13059.EP087938B02'),
    (13059, '13059.SP464352A03'),
    (13059, '13059.EP385387B01'),
    (13059, '13059.EP121013B01'),
    (13059, '13059.SP573849A04'),
    (13059, '13059.SP410796A02'),
    (13059, '13059.AP531397B04'),
    (13059, '13059.EP400448B04'),
    (13059, '13059.EP043583B01'),

    -- 6124
    (6124, '6124.GC.2atest.A'),
    (6124, '6124.T.Test.7.15.15B'),
    (6124, '6124.example.B.example'),
    (6124, '6124.Marine.Sediment.25.27cm.R2'),
    (6124, '6124.45272.1.swab.2'),
    (6124, '6124.T.Test.7.12.15A'),
    (6124, '6124.Test.4.3.2012.example'),
    (6124, '6124.Test.2.17.2014.example'),
    (6124, '6124.Marine.Sediment.20.22cm.R1'),
    (6124, '6124.Test.7.27.2014'),
    (6124, '6124.Test1.M1.B.1test.A'),
    (6124, '6124.Marine.Sediment.15.17cm.R1'),
    (6124, '6124.Test.4.4.2015.example'),
    (6124, '6124.T.Test.7.19.15A'),
    (6124, '6124.Test.9.28.2014.example'),
    (6124, '6124.Soil .Test.T2.2.Tube2'),
    (6124, '6124.Test.3.24.2015.example'),
    (6124, '6124.Test.12.28.2011.example'),
    (6124, '6124.Test1.T1.1test.A'),
    (6124, '6124.example.D.example'),
    (6124, '6124.T.Test.7.19.15B'),
    (6124, '6124.Test.10.27.2013.example'),
    (6124, '6124.45248.2.2'),
    (6124, '6124.Test.3.23.2014'),
    (6124, '6124.example.D.R2'),
    (6124, '6124.example.C.R4'),
    (6124, '6124.Soil.Test.T4.1.Tube4'),
    (6124, '6124.Test.2.8.2013'),
    (6124, '6124.Test.12.14.2013.example'),
    (6124, '6124.Test.6.29.2015.example'),
    (6124, '6124.45261.2.1'),
    (6124, '6124.A23'),
    (6124, '6124.Test.6.9.2013.example'),
    (6124, '6124.A30'),
    (6124, '6124.45327.7.2'),
    (6124, '6124.Test.8.10.2013.example'),
    (6124, '6124.T.Test.7.8.15A'),
    (6124, '6124.Marine.Sediment.0.2cm.R1'),
    (6124, '6124.1test.G.CNTL.A'),
    (6124, '6124.Soil.Test.T1.2.Tube1'),
    (6124, '6124.1test.M.CNTL.A'),
    (6124, '6124.Test.8.25.2014'),
    (6124, '6124.Test.1.14.2015'),
    (6124, '6124.T.Test.7.15.15B.example'),
    (6124, '6124.example.A.example'),
    (6124, '6124.Test.4.7.2013.example'),
    (6124, '6124.Test.1.6.2015.example'),
    (6124, '6124.Test.8.10.2013'),
    (6124, '6124.Test2.MT1.Sample.A'),
    (6124, '6124.Test.5.4.2014.example'),
    (6124, '6124.A21'),
    (6124, '6124.Test2.T1.01BH1.Y.A'),
    (6124, '6124.Test.7.26.2015.example'),
    (6124, '6124.45272.11.2'),
    (6124, '6124.Test1.T1.B.Sample.A'),
    (6124, '6124.Test.4.29.2013'),
    (6124, '6124.A27'),
    (6124, '6124.Test2.T1.B.A'),
    (6124, '6124.45326.1.swab.2'),
    (6124, '6124.example.B.R2'),
    (6124, '6124.Test.1.19.2014'),
    (6124, '6124.example.C.example'),
    (6124, '6124.Test.2.23.2015.example'),
    (6124, '6124.Test.11.6.2012.example'),
    (6124, '6124.Test.1.19.2014.example'),
    (6124, '6124.Test.8.22.2014.R1.example'),
    (6124, '6124.Test2.T1.A'),
    (6124, '6124.Marine.Sediment.10.12cm.R2'),
    (6124, '6124.Test.6.16.2014'),
    (6124, '6124.45316.8.1'),
    (6124, '6124.Test.8.22.2014.R2.example'),
    (6124, '6124.Test.1.26.2013'),
    (6124, '6124.Test.6.29.2015'),
    (6124, '6124.Test1.T1.A'),
    (6124, '6124.Marine.Sediment.5.7cm.R1'),
    (6124, '6124.Soil.Test.T4.2.Tube5'),
    (6124, '6124.Test.7.14.2013.example'),
    (6124, '6124.A31'),
    (6124, '6124.Test.11.16.2014'),
    (6124, '6124.Test.12.17.2014.example'),
    (6124, '6124.Test.2.25.2013.example'),
    (6124, '6124.Test2.MT1.1ex.Y.A'),
    (6124, '6124.Test.3.8.2015'),
    (6124, '6124.Soil.Test.T3.2.Tube3'),
    (6124, '6124.Test.9.3.2013.example'),
    (6124, '6124.Test.11.10.2013'),
    (6124, '6124.Test.3.24.2015'),
    (6124, '6124.45208.1.1'),
    (6124, '6124.Marine.Sediment.30.32cm.R3'),
    (6124, '6124.example.A.R3');
