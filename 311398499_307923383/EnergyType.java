
public enum EnergyType {
	EnergyWithoutEntropy, EnergyWithEntropy, EnergyForwarding;

	public static EnergyType phaseEnergy(String string) {
		switch (string.charAt(0)) {
		case '0':
			return EnergyWithoutEntropy;
		case '1':
			return EnergyWithEntropy;
		case '2':
			return EnergyForwarding;
		}
		return null;
	}
}
