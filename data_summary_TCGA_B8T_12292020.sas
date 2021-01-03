
***import b8T data to do analysis;
%inc "H:\Jerry\MACRO\char_to_num.sas";
proc import DBMS='CSV' datafile="H:\Rusty_B8T\data\TCGA_BLCA\data_b8t.csv" out=data_b8t replace;
     guessingrows=400;run;
%char_to_num(in_dsn=data_b8t,out_dsn=data_b8t2,var_list=death_days followup_days new_death);
proc sort data=data_b8t nodupkey; by shortid;run;

data data_b8t2;
     set data_b8t2;
	 new_death_m=new_death/30;
	 run;

proc format;
     value death 0='Censor'
	             1='Death';
	run;
ods rtf file="H:\Rusty_B8T\SAS\B8T_sample_data_summary_%sysfunc(today(),date9.).doc";
proc tabulate data=data_b8t2 missing f=8.;
     class b8t patient_gender death_event;
	 table patient_gender all, b8t*death_event all;
	 format death_event death.;
	 where new_death>=0;
	 title 'data summary of death and censor by gender and b8T';
	 run;

proc tabulate data=data_b8t2 missing f=8.1;
     class b8t patient_gender;
	 var new_death_m;
	 table patient_gender b8t all,new_death_m*(N*f=8. mean median min q1 q3 max);
	 where new_death>=0;
	 title 'data summary of follow-up time by gender and b8T';
	 run;


proc tabulate data=data_b8t2 missing f=8.1;
     class b8t patient_gender;
	 var new_death_m;
	 *table patient_gender b8t,new_death_m*(N*f=8. mean median min q1 q3 max);
	 table patient_gender*b8t all,new_death_m*(N*f=8. mean median min q1 q3 max);
	 where new_death>=0;
	 title 'data summary of follow-up time by gender and b8T';
	 run;
ods rtf close;
