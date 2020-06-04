import pandas as pd
import time

# EASY_PATH
# ============================================================================
def high_impact(df): 
    high_impact = ["splice_donor_variant","stop_gained","stop_gained, splice_region_variant",\
    "splice_acceptor_variant","stop_lost","frameshift_variant","frameshift_variant, feature_truncation",\
    "frameshift_variant, feature_elongation", "frameshift_variant, feature_truncation"]
    high_impact_variants = df[df["Consequence"].isin(high_impact)]
    return high_impact_variants 

def moderate_impact(df): 
    moderate_impact= ["missense_variant", "missense_variant, splice_region_variant",\
    "splice_region_variant, intron_variant","splice_region_variant, synonymous_variant","inframe_deletion","missense_variant,splice_region_variant"] 
    moderate_impact_variants = df[df["Consequence"].isin(moderate_impact)]
    return moderate_impact_variants
 

def in_silicoanalysis(df):
    moderate = moderate_impact(df)
    #PolyPhen2: Retrieving Scores 
    moderate.loc[:,"PolyPhen"] = moderate["PolyPhen"].str.extract('(\d*\.\d+|\d+)', expand=False).astype(float)
    #SIFT: Retrieving SIFT scores
    moderate.loc[:,"SIFT"] = moderate["SIFT"].str.extract('(\d*\.\d+|\d+)', expand=False).astype(float)                                       
    #Conditionals
    df_PolyPhen =moderate[moderate["PolyPhen"]>=0.40] 
    df_PP_SIFT = df_PolyPhen[df_PolyPhen["SIFT"] <= 0.05]
    return df_PP_SIFT
    

def clinical_data(df):
    ClinVar_Sig = ["pathogenic", "unknown","uncertain_significance", "other:pathogenic"]
    high = high_impact(df)
    moderate =in_silicoanalysis(df)
    if high.empty and not moderate.empty:
        moderate_annotated = moderate[moderate["CLIN_SIG"].isin(ClinVar_Sig)]
        if moderate_annotated.empty:
            return moderate
        else:
            return moderate_annotated
    if moderate.empty and not high.empty:
            high_annotated = high[high["CLIN_SIG"].isin(ClinVar_Sig)]
            if high_annotated.empty:
                return high
            else:
                return high_annotated  
    if high.empty and moderate.empty:
        return "Your file seems not to have potential candidate variants. You may need to check it manually"
    else:
        merge = pd.concat([high, moderate])
        merge_annotated = merge[merge["CLIN_SIG"].isin(ClinVar_Sig)]
        if merge_annotated.empty:
            return merge
        else: 
            return merge_annotated
        

def easy_path():
    df = pd.read_excel(input("Please enter the Pre-annotated file name: "))
    start_time = time.time()
    output = ["Gene","HGVSc","Consequence","Existing_variation", "CLIN_SIG", "PHENOTYPES"]
    results = clinical_data(df)
    end_time = time.time()
    run_time = "Easy Path took",end_time-start_time,"seconds in running your file"
    if isinstance(results, pd.DataFrame):
        results_output = results.loc[:,results.columns.isin(output)]
        return results_output, run_time
    else: 
        return results, run_time
    

#=============================================================================
#Testing
print(easy_path())
