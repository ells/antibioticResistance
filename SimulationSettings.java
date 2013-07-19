package antibioticResistance;

/**
 * Created with IntelliJ IDEA.
 * User: ellscampbell
 * Date: 7/16/13
 * Time: 4:36 PM
 * To change this template use File | Settings | File Templates.
 */
public class SimulationSettings {

    public static final int uninfected = 0;
    public static final int susceptible = 1;
    public static final int resistantA = 2;
    public static final int resistantB = 3;
    public static final int multiResistant = 4;

    public static final int noTreatment = 0;
    public static final int singleTreatment = 1;
    public static final int separateCocktail = 2;
    public static final int combinedCocktail = 3;

    private int currentTreatment = noTreatment;

    private static SimulationSettings ourInstance = new SimulationSettings();
    private int period = 50;

    private double cocktailDoseCoefficient = 0.5;
    private int singleDosageA = 240;
    private int singleDosageB = 240;

    private double cocktailDosageA = singleDosageA * cocktailDoseCoefficient;
    private double cocktailDosageB = singleDosageB * cocktailDoseCoefficient;
    
    private double cost = 0.1;
    private double tradeoff = 1;

    private double rateOfInfection = 1;
    private double immuneClearance = 0.25;
    private double mutationRate = 0.000001;
    private double rateOfSuperInfection = 0.25;

    private double maximumBacterialGrowthRate = 0.25;
    private double minimumBacterialGrowthRate_WithMaxDrug = -0.25;

    private int baseResistantMIC = 240;
    private double baseSusceptibleMIC = baseResistantMIC * 0.000001;

    private double micSusceptible_drugA = baseSusceptibleMIC;
    private double micSusceptible_drugB = baseSusceptibleMIC;
    
    private double micResistantA_drugA = baseResistantMIC;
    private double micResistantA_drugB = baseSusceptibleMIC;

    private double micResistantB_drugA = baseSusceptibleMIC;
    private double micResistantB_drugB = baseResistantMIC;

    private double micMultiResistant_drugA;
    private double micMultiResistant_drugB;
    
    private double healthyPatients = 0.8;
    private double susceptible_InfectedPatients = 0.2;
    private double resistantA_InfectedPatients = 0.0;
    private double resistantB_InfectedPatients = 0.0;
    private double multiResistant_InfectedPatients = 0.0;

    public static SimulationSettings getInstance() {
        return ourInstance;
    }

    private SimulationSettings() {
    }

    public void setSingleDosageA(int dosageA) {
        this.singleDosageA = dosageA;
    }

    public void setMutationRate(double mu) {
        this.mutationRate = mu;
    }

    public void setSingleDosageB(int dosageB) {
        this.singleDosageB = dosageB;
    }

    public int getSingleDosageA() {
        return singleDosageA;
    }

    public int getSingleDosageB() {
        return singleDosageB;
    }

    public int getCurrentTreatment() {
        return currentTreatment;
    }

    public double getBaseSusceptibleMIC() {
        return baseSusceptibleMIC;
    }

    public int getBaseResistantMIC() {
        return baseResistantMIC;
    }

    public double getCocktailDosageA() {
        return cocktailDosageA;
    }

    public double getCocktailDosageB() {
        return cocktailDosageB;
    }

    public int getPeriod() {
        return period;
    }

    public double getImmuneClearance() {
        return immuneClearance;
    }

    public double getMutationRate() {
        return mutationRate;
    }

    public double getRateOfInfection() {
        return rateOfInfection;

    }

    public double getCost() {
        return cost;
    }

    public double getTradeoff() {
        return tradeoff;
    }

    public double getRateOfSuperInfection() {
        return rateOfSuperInfection;
    }

    public double getCocktailDoseCoefficient() {
        return cocktailDoseCoefficient;
    }

    public void setTradeoff(double tradeoff) {
        this.tradeoff = tradeoff;
        this.micMultiResistant_drugA = baseResistantMIC * cocktailDoseCoefficient * (this.tradeoff);
        this.micMultiResistant_drugB = baseResistantMIC * cocktailDoseCoefficient * (this.tradeoff);
    }


    public double getMinimumBacterialGrowthRate_WithMaxDrug() {
        return minimumBacterialGrowthRate_WithMaxDrug;
    }

    public double getMaximumBacterialGrowthRate() {
        return maximumBacterialGrowthRate;
    }

    public double getMicSusceptible_drugA() {
        return micSusceptible_drugA;
    }

    public double getMicSusceptible_drugB() {
        return micSusceptible_drugB;
    }

    public double getMicResistantA_drugA() {
        return micResistantA_drugA;
    }

    public double getMicResistantA_drugB() {
        return micResistantA_drugB;
    }

    public double getMicResistantB_drugA() {
        return micResistantB_drugA;
    }

    public double getMicResistantB_drugB() {
        return micResistantB_drugB;
    }

    public double getMicMultiResistant_drugA() {
        return micMultiResistant_drugA;
    }

    public double getMicMultiResistant_drugB() {
        return micMultiResistant_drugB;
    }


    public double getHealthyPatients() {
        return healthyPatients;
    }

    public double getSusceptible_InfectedPatients() {
        return susceptible_InfectedPatients;
    }

    public double getResistantA_InfectedPatients() {
        return resistantA_InfectedPatients;
    }

    public double getResistantB_InfectedPatients() {
        return resistantB_InfectedPatients;
    }

    public double getMultiResistant_InfectedPatients() {
        return multiResistant_InfectedPatients;
    }




}
