% batch job

job=batch('periodic_Dlt_SphHMC');

wait(job);
diary(job,'periodic_Dlt_SphHMC_diary');
load(job);

delete(job);
clear job;