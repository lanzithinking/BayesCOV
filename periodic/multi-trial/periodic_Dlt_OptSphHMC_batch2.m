% batch job

job=batch('periodic_Dlt_OptSphHMC2');

wait(job);
diary(job,'periodic_Dlt_OptSphHMC_diary2');
load(job);

delete(job);
clear job;