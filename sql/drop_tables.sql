-- $Id: drop_tables.sql,v 1.1 2013-01-29 12:12:51 ak4 Exp $

-- This SQL file drops all the tables
-- created bf "table.sql", if they exist.

DROP TABLE IF EXISTS job;
DROP TABLE IF EXISTS job_status;
DROP TABLE IF EXISTS rule_goal;
DROP TABLE IF EXISTS rule_conditions;
DROP TABLE IF EXISTS input_id_analysis;
DROP TABLE IF EXISTS input_id_type_analysis;
