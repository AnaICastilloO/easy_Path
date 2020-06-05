
#Modules used for Easy_Path and Performance Test 
import pandas as pd
import time
import matplotlib.pyplot as plt


# EASY_PATH
# ============================================================================

# High_Impact Function 
def high_impact(df): 
    """
    Args:
        dataframe
    Returns:
        new dataframe with variants associated to a hight impact
    """
    #High Impact list: 
    high_impact = ["splice_donor_variant","stop_gained","stop_gained, splice_region_variant",\
    "splice_acceptor_variant","stop_lost","frameshift_variant","frameshift_variant, feature_truncation",\
    "frameshift_variant, feature_elongation", "frameshift_variant, feature_truncation"]
    #Selecting variants with high consequence
    high_impact_variants = df[df["Consequence"].isin(high_impact)]
    return high_impact_variants 

#Moderate Function
def moderate_impact(df): 
    """
    Args:
        dataframe
    Returns:
        new dataframe with variants associated to a moderate impact
    """
    #Moderate list:
    moderate_impact= ["missense_variant", "missense_variant, splice_region_variant",\
    "splice_region_variant, intron_variant","splice_region_variant, synonymous_variant",\
        "inframe_deletion","missense_variant,splice_region_variant"] 
    #Selecting variants with moderate consequence
    moderate_impact_variants = df[df["Consequence"].isin(moderate_impact)]
    return moderate_impact_variants
 
# In_silico analysis Function
def in_silicoanalysis(df):
     """
    Args:
        dataframe
    Returns:
        data frame with variants passing scores with PolyPhen and Sift accepted thresholds
    """ 
     moderate = moderate_impact(df)
    #PolyPhen2: Retrieving Scores 
     moderate.loc[:,"PolyPhen"] = moderate["PolyPhen"].str.extract('(\d*\.\d+|\d+)', \
     expand=False).astype(float)
    #SIFT: Retrieving SIFT scores
     moderate.loc[:,"SIFT"] = moderate["SIFT"].str.extract('(\d*\.\d+|\d+)', \
     expand=False).astype(float)                                       
    #Conditionals
     df_PolyPhen =moderate[moderate["PolyPhen"]>=0.40] 
     df_PP_SIFT = df_PolyPhen[df_PolyPhen["SIFT"] < 0.05]
     return df_PP_SIFT
    
# Clinical_Data Function
def clinical_data(df):
     """
    Args:
        dataframe
    Returns:
        Depending of high and moderate function outputs, this functions returns a 
        combination of dataframes from high or moderate impact, highlightin those with clinical
        annotation. In case of absense, this function will return a string that not potential variants 
        were identified. 
    """ 
    #List with Clinical Significance. Relevant to EasyPath Program
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
        return "Your file seems not to have potential candidate variants. \
            You may need to check it manually"
     else:
        merge = pd.concat([high, moderate])
        merge_annotated = merge[merge["CLIN_SIG"].isin(ClinVar_Sig)]
        if merge_annotated.empty:
            return merge
        else: 
            return merge_annotated
        
# Easy_Path Function (EASY_PATH)
def easy_path():
    """
    Args:
        no arguments
    Returns:
        This function returns the final dataframe encompassing relevant 
        information for the clinician and the run time Easy Path took in running 
        the pre-annotated file 
    """ 
     
    #Runnning times and calling clinical_data function
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




#=============================================================================
#PERFORMANCE TEST 
Files = ["T141", "T03", "T59", "T26", "T33", "T225", "T08", "T47", "T38", "T48", "T06", "T52", "T134"]
time = [174, 180, 210, 219, 222, 225, 228, 228, 246, 258, 270, 312, 324]
# plotting the line 1 points 
plt.plot(Files, time, label = "manual", marker="D", color="green")
# line 2 points
Files = ["T141", "T03", "T59", "T26", "T33", "T225", "T08", "T47", "T38", "T48", "T06", "T52", "T134"]
time = [0.005, 0.014, 0.004, 0.004, 0.004, 0.014, 0.014, 0.003, 0.014, 0.015, 0.016, 0.014, 0.016]
# plotting the line 2 points 
plt.plot(Files, time, label = "python", marker="*", color="purple")
plt.xlabel('ExcelFiles')
# Set the y axis label of the current axis.
plt.ylabel('Time in seconds')
# Set a title of the current axes.
# to increase the quality of image
plt.savefig('Time difference between manual and python methods', dpi=300)
plt.title('Time differences between manual and python methods')
# show a legend on the plot
plt.legend()
# Display a figure.
plt.show()
