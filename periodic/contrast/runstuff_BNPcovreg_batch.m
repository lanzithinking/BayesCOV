% batch job

job=batch('runstuff_BNPcovreg');

wait(job);
diary(job,'runstuff_BNPcovreg_diary');
load(job);

delete(job);
clear job;