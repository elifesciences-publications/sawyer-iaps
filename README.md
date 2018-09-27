# Code and Data for "Alcoholism Gender Differences in Brain Responsivity to Emotional Stimuli"

This repository contains the code and preprocessed data used to generate tables and figures for
Sawyer KS, Maleki N, Urban T, Marinkovic K, Karson S, Ruiz SM, Harris GJ, 
"Alcoholism Gender Differences in Brain Responsivity to Emotional Stimuli"
(https://www.biorxiv.org/content/early/2018/12/27/428565)

The code ("R_code/Sawyer-IAPS.R") was written to work with R version 3.5.0.

The preprocessed data ("R_data/Sawyer-IAPS.csv") are in comma separated value format. The columns represent various features, including demographic attributes, participant responses to stimuli, and functional magnetic resonance imaging (fMRI) contrast values, as described in the methods of the associated manuscript. For the participant responses measures, an example that describes how the variables were coded follows.

pics_percent_pos_bad: percent of responces to positive ("happy") pictures rated bad. 

The following abbreviations were used for this coding system.

pics: pictures

percent: percent pictures rated

RT: reaction time

pos: positive or "happy" pictures

neg: negative or "aversive" pictures

neutral: neutral pictures

erotic: erotic pictures

gross: "gruesome" pictures

bad: bad ratings

good: good ratings

neut: neutral ratings



Further variables are as follows.

ID: participant number

Age: rounded to nearest year

VIQ: Verbal Intelligence Quotient from Wechsler Adult Intelligence Scale

PIQ: Performance Intelligence Quotient from Wechsler Adult Intelligence Scale

FSIQ: Full Scale Intelligence Quotient from Wechsler Adult Intelligence Scale

WMS_IMI: Immediate Memory Index from Wechsler Memory Scale

WMS_DMI: Delayed (General) Memory Index from Wechsler Memory Scale

WMS_WMI: Working Memory Index from Wechsler Memory Scale

HRSD: Hamilton Rating Scale for Depression

DHD: Duration of Heavy Drinking (21 or more drinks per week), in years

DD: Daily Drinks, in ounces of ethanol per day

LOS: Length of Sobriety, in years

Ten: Profile of Mood States Tension

Dep: Profile of Mood States Depression

Ang: Profile of Mood States Anger

Vig: Profile of Mood States Vigor

Fat: Profile of Mood States Fatigue

Con: Profile of Mood States Confusion

A: Multiple Adjective Affect Check List Anxiety

D: Multiple Adjective Affect Check List Depression

H: Multiple Adjective Affect Check List Hostility

PA: Multiple Adjective Affect Check List Positive Affect

SS: Multiple Adjective Affect Check List Sensation Seeking

Dys: Multiple Adjective Affect Check List Dysphoria

PASS: Multiple Adjective Affect Check List Positive Affect Sensation Seekingc

behavioral_responses_included: '1' indicates participants for which behavioral responses were available

Limbic Structures: Erotic vs neutral contrast effect size values from volumetric cluster

caudalmiddlefrontal1: Negative vs neutral contrast effect size values from surface cluster 1 with peak vertex in caudal middle frontal cortex

caudalmiddlefrontal2: Negative vs neutral contrast effect size values from surface cluster 2 with peak vertex in caudal middle Frontal cortex

superiorfrontal: Negative vs neutral contrast effect size values from surface cluster 2 with peak vertex in caudal superior frontal cortex

inferiorparietal: Negative vs neutral contrast effect size values from surface cluster 2 with peak vertex in inferior parietal cortex
