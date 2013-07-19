package antibioticResistance;

import java.util.ArrayList;

/**
 * Created with IntelliJ IDEA.
 * User: ellscampbell
 * Date: 7/16/13
 * Time: 4:36 PM
 * To change this template use File | Settings | File Templates.
 */
public class staticSimulation {

    // variable reference and storage
    private static final int micA = 0;
    private static final int micB = 1;

    private static final int uninfected = 0;
    private static final int susceptible = 1;
    private static final int resistantA = 2;
    private static final int resistantB = 3;
    private static final int multiResistant = 4;

    private static final int noTreatment = 0;
    private static final int singleTreatmentA = 1;
    private static final int singleTreatmentB = 2;
    private static final int separateCocktail = 3;
    private static final int combinedCocktail = 4;

    double[] pharmacoGrowth = new double[5];
    double[] patients = new double[5];
    private int deltaTime = 1;
    private double PsiMax, PsiMin, mu, immuneClearance, rateOfInfection, rateOfSuperInfection, cost;
    private int dosageA, dosageB;


    public static void main(String[] args) {
        staticSimulation sim = new staticSimulation();
        sim.run();
    }

    public void run() {
        for (int timeStep = 0; timeStep < 100000; timeStep++) {
            fourthOrderApproximation();
        }
        System.out.println((1.0 - (SimulationSettings.getInstance().getTradeoff()/2.0)) + "\t" + patients[uninfected] + "\t" + patients[susceptible] + "\t" + patients[resistantA] + "\t" + patients[resistantB] + "\t" + patients[multiResistant]);

    }

    public void initConstants() {
        // bergstrom constants
        cost = SimulationSettings.getInstance().getCost();
        immuneClearance = SimulationSettings.getInstance().getImmuneClearance();
        rateOfInfection = SimulationSettings.getInstance().getRateOfInfection();
        rateOfSuperInfection = SimulationSettings.getInstance().getRateOfSuperInfection();

        // pharmacodynamic constants
        PsiMax = SimulationSettings.getInstance().getMaximumBacterialGrowthRate();
        PsiMin = SimulationSettings.getInstance().getMinimumBacterialGrowthRate_WithMaxDrug();

        // our constants
        mu = SimulationSettings.getInstance().getMutationRate();
        dosageA = SimulationSettings.getInstance().getSingleDosageA();
        dosageB = SimulationSettings.getInstance().getSingleDosageB();
    }



    private void fourthOrderApproximation() {
        ArrayList<double[]> fourDeltas = new ArrayList<double[]>() ;

        for (int i = 0; i < 3; i++) {
            double[] deltaPatients = new double[5];
            deltaPatients[uninfected] = deltaUninfected();
            deltaPatients[susceptible] = deltaSusceptible();
            deltaPatients[resistantA] = deltaResistantA();
            deltaPatients[resistantB] = deltaResistantB();
            deltaPatients[multiResistant] = deltaMultiResistant();

            patients[uninfected]     += 0.5 * deltaTime * deltaPatients[uninfected];
            patients[susceptible]    += 0.5 * deltaTime * deltaPatients[susceptible];
            patients[resistantA]     += 0.5 * deltaTime * deltaPatients[resistantA];
            patients[resistantB]     += 0.5 * deltaTime * deltaPatients[resistantB];
            patients[multiResistant] += 0.5 * deltaTime * deltaPatients[multiResistant];

            fourDeltas.add(deltaPatients);
        }

        double[] deltaPatients = new double[5];

        deltaPatients[uninfected]     = deltaTime * deltaUninfected();
        deltaPatients[susceptible]    = deltaTime * deltaSusceptible();
        deltaPatients[resistantA]     = deltaTime * deltaResistantA();
        deltaPatients[resistantB]     = deltaTime * deltaResistantB();
        deltaPatients[multiResistant] = deltaTime * deltaMultiResistant();

        fourDeltas.add(deltaPatients);



        patients[uninfected]     += (fourDeltas.get(0)[uninfected] + (2.0 * fourDeltas.get(1)[uninfected]) + (2.0 * fourDeltas.get(2)[uninfected]) + fourDeltas.get(3)[uninfected]) * ((1.0/6.0) * deltaTime);
        patients[susceptible]    += (fourDeltas.get(0)[susceptible] + (2.0 * fourDeltas.get(1)[susceptible]) + (2.0 * fourDeltas.get(2)[susceptible]) + fourDeltas.get(3)[susceptible]) * ((1.0/6.0) * deltaTime);
        patients[resistantA]     += (fourDeltas.get(0)[resistantA] + (2.0 * fourDeltas.get(1)[resistantA]) + (2.0 * fourDeltas.get(2)[resistantA]) + fourDeltas.get(3)[resistantA]) * ((1.0/6.0) * deltaTime);
        patients[resistantB]     += (fourDeltas.get(0)[resistantB] + (2.0 * fourDeltas.get(1)[resistantB]) + (2.0 * fourDeltas.get(2)[resistantB]) + fourDeltas.get(3)[resistantB]) * ((1.0/6.0) * deltaTime);
        patients[multiResistant] += (fourDeltas.get(0)[multiResistant] + (2.0 * fourDeltas.get(1)[multiResistant]) + (2.0 * fourDeltas.get(2)[multiResistant]) + fourDeltas.get(3)[multiResistant]) * ((1.0/6.0) * deltaTime);



    }

    public void initFrequencies() {
        this.patients[uninfected] = SimulationSettings.getInstance().getHealthyPatients();
        this.patients[susceptible] = SimulationSettings.getInstance().getSusceptible_InfectedPatients();
        this.patients[resistantA] = SimulationSettings.getInstance().getResistantA_InfectedPatients();
        this.patients[resistantB] = SimulationSettings.getInstance().getResistantB_InfectedPatients();
        this.patients[multiResistant] = SimulationSettings.getInstance().getMultiResistant_InfectedPatients();
    }

    public void initCombinedGrowth() {
        this.pharmacoGrowth[uninfected]     = 0;
        this.pharmacoGrowth[susceptible]    = calculateNetGrowthCOMBINED(susceptible);
        this.pharmacoGrowth[resistantA]     = calculateNetGrowthCOMBINED(resistantA);
        this.pharmacoGrowth[resistantB]     = calculateNetGrowthCOMBINED(resistantB);
        this.pharmacoGrowth[multiResistant] = calculateNetGrowthCOMBINED(multiResistant);

    }

    public void initSeparateGrowth() {
        this.pharmacoGrowth[uninfected]     = 0;
        this.pharmacoGrowth[susceptible]    = calculateNetGrowthSEPARATE(susceptible);
        this.pharmacoGrowth[resistantA]     = calculateNetGrowthSEPARATE(resistantA);
        this.pharmacoGrowth[resistantB]     = calculateNetGrowthSEPARATE(resistantB);
        this.pharmacoGrowth[multiResistant] = calculateNetGrowthSEPARATE(multiResistant);
    }

    public void initSingleAGrowth() {
        this.pharmacoGrowth[uninfected]     = 0;
        this.pharmacoGrowth[susceptible]    = calculateNetGrowth_SingleA(susceptible);
        this.pharmacoGrowth[resistantA]     = calculateNetGrowth_SingleA(resistantA);
        this.pharmacoGrowth[resistantB]     = calculateNetGrowth_SingleA(resistantB);
        this.pharmacoGrowth[multiResistant] = calculateNetGrowth_SingleA(multiResistant);
    }

    public void initSingleBGrowth() {
        this.pharmacoGrowth[uninfected]     = 0;
        this.pharmacoGrowth[susceptible]    = calculateNetGrowth_SingleB(susceptible);
        this.pharmacoGrowth[resistantA]     = calculateNetGrowth_SingleB(resistantA);
        this.pharmacoGrowth[resistantB]     = calculateNetGrowth_SingleB(resistantB);
        this.pharmacoGrowth[multiResistant] = calculateNetGrowth_SingleB(multiResistant);
    }

    public void initNoTreatGrowth() {
        this.pharmacoGrowth[uninfected]     = 0;
        this.pharmacoGrowth[susceptible]    = calculateNetGrowth_NoTreat();
        this.pharmacoGrowth[resistantA]     = calculateNetGrowth_NoTreat();
        this.pharmacoGrowth[resistantB]     = calculateNetGrowth_NoTreat();
        this.pharmacoGrowth[multiResistant] = calculateNetGrowth_NoTreat();
    }

    private double calculateNetGrowthCOMBINED(int strain) {
         double[] MICs = selectMICs(strain);
         double micA = MICs[this.micA];
         double micB = MICs[this.micB];

         double netGrowth, combinedEffect;
         combinedEffect = (PsiMax - PsiMin) * ((dosageA/micA) + (dosageB/micB)) / (((dosageA/micA) + (dosageB/micB)) - (PsiMin/PsiMax));
         netGrowth = PsiMax - combinedEffect;
         return netGrowth;
    }


    private double calculateNetGrowth_NoTreat() {
        return PsiMax;
    }

    private double calculateNetGrowth_SingleA(int strain) {
        double[] MICs = selectMICs(strain);
        double micA = MICs[this.micA];
        double micB = MICs[this.micB];
        double netGrowth = 0;

        double effectA = ((PsiMax - PsiMin) * ((dosageA/ micA ))) / (((dosageA/ micA )) - (PsiMin / PsiMax));

        netGrowth = PsiMax - effectA;

        return netGrowth;
    }

    private double calculateNetGrowth_SingleB(int strain) {
        double[] MICs = selectMICs(strain);
        double micA = MICs[this.micA];
        double micB = MICs[this.micB];
        double netGrowth = 0;

        double effectB = ((PsiMax - PsiMin) * ((dosageB/ micB ))) / (((dosageB/ micB )) - (PsiMin / PsiMax));
        netGrowth = PsiMax - effectB;
        return netGrowth;
    }

    private double[] selectMICs(int strain) {
        double[] MICs = new double[2];
        double micA = 0;
        double micB = 0;

        if (strain == susceptible) {
            micA = SimulationSettings.getInstance().getMicSusceptible_drugA();
            micB = SimulationSettings.getInstance().getMicSusceptible_drugB();
        }

        if (strain == resistantA) {
            micA = SimulationSettings.getInstance().getMicResistantA_drugA();
            micB = SimulationSettings.getInstance().getMicResistantA_drugB();
        }

        if (strain == resistantB) {
            micA = SimulationSettings.getInstance().getMicResistantB_drugA();
            micB = SimulationSettings.getInstance().getMicResistantB_drugB();
        }

        if (strain == multiResistant) {
            micA = SimulationSettings.getInstance().getMicMultiResistant_drugA();
            micB = SimulationSettings.getInstance().getMicMultiResistant_drugB();
        }

        MICs[this.micA] = micA;
        MICs[this.micB] = micB;
        return MICs;
    }

    private double calculateNetGrowthSEPARATE(int strain) {
        double[] MICs = selectMICs(strain);
        double micA = MICs[this.micA];
        double micB = MICs[this.micB];

        double netGrowth, effectA, effectB;

        effectA = ((PsiMax - PsiMin) * ((dosageA/ micA ))) / (((dosageA/ micA )) - (PsiMin / PsiMax));
        effectB = ((PsiMax - PsiMin) * ((dosageB/ micB ))) / (((dosageB/ micB )) - (PsiMin / PsiMax));

        netGrowth = PsiMax - effectA - effectB;

        return netGrowth;
    }


    private double deltaSusceptible() {
        double deltaS = 0;
        double mutation, clearance, infection, superInfection;


        mutation =  (mu * patients[resistantA]) + (mu * patients[resistantB]) + (mu * patients[multiResistant]) - (3 * mu * patients[susceptible]);
        clearance =  (pharmacoGrowth[susceptible] - immuneClearance) * patients[susceptible];
        infection =  patients[susceptible] * patients[uninfected] * rateOfInfection;
        superInfection = (patients[susceptible] * rateOfSuperInfection * rateOfInfection * cost) * (patients[resistantA] + patients[resistantB] + patients[multiResistant]);

        deltaS = mutation + clearance + infection + superInfection ;

        return deltaS;
    }

    private double deltaResistantA() {
        double deltaResistantA = 0;
        double mutation, clearance, infection, superInfection;

        mutation = (mu * patients[susceptible]) + (mu * patients[resistantB]) + (mu * patients[multiResistant]) - (3 * mu * patients[resistantA]);
        clearance = (pharmacoGrowth[resistantA] - immuneClearance) * patients[resistantA];
        infection =  (patients[resistantA] * (1.0-cost)) * patients[uninfected] * rateOfInfection;
        superInfection = patients[resistantA] * rateOfSuperInfection * rateOfInfection * cost * patients[susceptible];

        deltaResistantA = mutation + clearance + infection - superInfection;

        return deltaResistantA;
    }

    private double deltaResistantB() {
        double deltaResistantB = 0;
        double mutation, clearance, infection, superInfection;


        mutation = (mu * patients[susceptible]) + (mu * patients[resistantA]) + (mu * patients[multiResistant]) - (3 * mu * patients[resistantB]);
        clearance = (pharmacoGrowth[resistantB] - immuneClearance) * patients[resistantB];
        infection =  (patients[resistantB] * (1.0-cost)) * patients[uninfected] * rateOfInfection;
        superInfection = patients[resistantB] * rateOfSuperInfection * rateOfInfection * cost * patients[susceptible];

        deltaResistantB = mutation + clearance + infection - superInfection;

        return deltaResistantB;
    }

    private double deltaMultiResistant() {
        double deltaMultiResistant = 0;
        double mutation, clearance, infection, superInfection;


        mutation = (mu * patients[susceptible]) + (mu * patients[resistantA]) + (mu * patients[resistantB]) - (3 * mu * patients[multiResistant]);
        clearance = (pharmacoGrowth[multiResistant] - immuneClearance) * patients[multiResistant];
        infection =  (patients[multiResistant] * (1.0-cost)) * patients[uninfected] * rateOfInfection;
        superInfection = patients[multiResistant] * rateOfSuperInfection * rateOfInfection * cost * patients[susceptible];

        deltaMultiResistant = mutation + clearance + infection - superInfection;

        return deltaMultiResistant;
    }

    private double deltaUninfected() {
        double deltaUninfected = 0;

        deltaUninfected = deltaUninfected - (pharmacoGrowth[susceptible] - immuneClearance) * patients[susceptible];
        deltaUninfected = deltaUninfected - (pharmacoGrowth[resistantA] - immuneClearance) * patients[resistantA];
        deltaUninfected = deltaUninfected - (pharmacoGrowth[resistantB] - immuneClearance) * patients[resistantB];
        deltaUninfected = deltaUninfected - (pharmacoGrowth[multiResistant] - immuneClearance) * patients[multiResistant];

        deltaUninfected = deltaUninfected - (patients[susceptible] * rateOfInfection * patients[uninfected]);
        deltaUninfected = deltaUninfected - (patients[resistantA] * rateOfInfection * (1.0 - cost) * patients[uninfected]) ;
        deltaUninfected = deltaUninfected - (patients[resistantB] * rateOfInfection * (1.0 - cost) * patients[uninfected]) ;
        deltaUninfected = deltaUninfected - (patients[multiResistant] * rateOfInfection * (1.0 - cost) * patients[uninfected]);

        return deltaUninfected;

    }


}
