package antibioticResistance;

/**
 * Created with IntelliJ IDEA.
 * User: ellscampbell
 * Date: 7/16/13
 * Time: 4:35 PM
 * To change this template use File | Settings | File Templates.
 */
public class SimulationManager {

    private static final int noTreatmentControl = 0;
    private static final int singleTreatmentA = 1;
    private static final int singleTreatmentB = 2;
    private static final int separateCocktail = 3;
    private static final int combinedCocktail = 4;
    private static final int cycle = 5;
    private static final int mix = 6;
    public static final int noMutationControl = 7;


    public static void main(String[] args) {
        SimulationManager sm = new SimulationManager();
        sm.exploreTradeoff(singleTreatmentB);
    }

    private void exploreTradeoff(int treatment) {

        if (treatment == mix) System.out.println("tradeoff\tuninfected\tsusceptible_drugA\tsusceptible_drugB\tresistantA_drugA\tresistantA_drugB\tresistantB_drugA\tresistantB_drugB\tmultiResist_drugA\tmultiResist_drugB");
        else System.out.println("tradeoff\tuninfected\tsusceptible\tresistantA\tresistantB\tmultiResistant");


        for (int i = 1; i < 100; i++) {
            double tradeoff = (i) * 0.02;
            SimulationSettings.getInstance().setTradeoff(tradeoff);
            runSpecifiedTreatment(treatment);
        }
    }

    private void runSpecifiedTreatment(int treatment) {
        if (treatment == noTreatmentControl) runNoTreatControl();
        if (treatment == singleTreatmentA)   runSingleA();
        if (treatment == singleTreatmentB)   runSingleB();
        if (treatment == separateCocktail)   runSeparateCocktail();
        if (treatment == combinedCocktail)   runCombinedCocktail();
        if (treatment == cycle)              runCycle();
        if (treatment == mix)                runMix();
        if (treatment == noMutationControl)  runNoMutationControl();
    }

    private void runCycle() {
        cycleSimulation cycle = new cycleSimulation();
        cycle.initFrequencies();
        cycle.initConstants();
        cycle.run();
    }

    private void runMix() {
        mixedSimulation mix = new mixedSimulation();
        mix.initFrequencies();
        mix.initConstants();
        mix.run();
    }

    private void runSeparateCocktail() {
        SimulationSettings.getInstance().setSingleDosageA(120);
        SimulationSettings.getInstance().setSingleDosageB(120);
        staticSimulation staticSim = new staticSimulation();
        staticSim.initFrequencies();
        staticSim.initConstants();
        staticSim.initSeparateGrowth();
        staticSim.run();
    }

    private void runCombinedCocktail() {
        SimulationSettings.getInstance().setSingleDosageA(120);
        SimulationSettings.getInstance().setSingleDosageB(120);
        staticSimulation staticSim = new staticSimulation();
        staticSim.initFrequencies();
        staticSim.initConstants();
        staticSim.initCombinedGrowth();
        staticSim.run();
    }

    private void runSingleA() {
        SimulationSettings.getInstance().setSingleDosageA(240);
        SimulationSettings.getInstance().setSingleDosageB(0);
        staticSimulation staticSim = new staticSimulation();
        staticSim.initFrequencies();
        staticSim.initConstants();
        staticSim.initSingleAGrowth();
        staticSim.run();
    }

    private void runNoMutationControl() {
        SimulationSettings.getInstance().setSingleDosageA(240);
        SimulationSettings.getInstance().setSingleDosageB(0);
        SimulationSettings.getInstance().setMutationRate(0.0);
        staticSimulation staticSim = new staticSimulation();
        staticSim.initFrequencies();
        staticSim.initConstants();
        staticSim.initSingleAGrowth();
        staticSim.run();
    }

    private void runSingleB() {
        SimulationSettings.getInstance().setSingleDosageA(0);
        SimulationSettings.getInstance().setSingleDosageB(240);
        staticSimulation staticSim = new staticSimulation();
        staticSim.initFrequencies();
        staticSim.initConstants();
        staticSim.initSingleBGrowth();
        staticSim.run();
    }

    private void runNoTreatControl() {
        SimulationSettings.getInstance().setSingleDosageA(0);
        SimulationSettings.getInstance().setSingleDosageB(0);
        staticSimulation staticSim = new staticSimulation();
        staticSim.initFrequencies();
        staticSim.initConstants();
        staticSim.initNoTreatGrowth();
        staticSim.run();
    }

}
