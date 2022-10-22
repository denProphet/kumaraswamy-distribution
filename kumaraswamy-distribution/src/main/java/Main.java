import distribution.Interval;
import distribution.KumaraswamyDistribution;
import distribution.exceptions.WrongIntervalException;

public class Main {
    public static void main(String[] args) {
        Interval interval = new Interval(0.5,0.5);
        KumaraswamyDistribution kumaraswamyDistribution = new KumaraswamyDistribution(interval);
        try {
            System.out.println(kumaraswamyDistribution.ppf(0.5));
        }catch (WrongIntervalException e){
            e.printStackTrace();
        }
        System.out.println(kumaraswamyDistribution.dispersion());
    }
}
