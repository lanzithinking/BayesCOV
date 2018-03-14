% batch job

job=batch('periodic_Dlt_OptSphHMC');

wait(job);
diary(job,'periodic_Dlt_OptSphHMC_diary');
load(job);

delete(job);
clear job;