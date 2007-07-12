
# this patch needs to be applied to reference databases if they used an ensembl-pipeline 
# checkout before January 2007 . 'condition' has become a reseved word in mysql 
# ( changed in january 2007 ) 
 
ALTER TABLE rule_conditions CHAGE condition rule_condition VARCHAR(40)  ; 


