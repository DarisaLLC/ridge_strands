
/**
 * The Class Correction.
 */
class Correction {

	/** The w est. */
	float w_est; /* Total line width extracted from the image */

	/** The r est. */
	float r_est; /* Gradient ratio extracted from the image */

	/** The w. */
	float w; /* True line width */

	/** The h. */
	float h; /* True asymmetry */

	/** The correction. */
	float correction; /* Line position correction */

	/** The w strong. */
	float w_strong; /* True width on the side with the stronger gradient */

	/** The w weak. */
	float w_weak; /* True width on the side with the weaker gradient */

	/** The is valid. */
	bool is_valid; /* Is this table entry valid? */

	/**
	 * Instantiates a new correction.
	 *
	 * @param w_est
	 *            the w est
	 * @param r_rest
	 *            the r rest
	 * @param w
	 *            the w
	 * @param h
	 *            the h
	 * @param correction
	 *            the correction
	 * @param w_strong
	 *            the w strong
	 * @param w_weak
	 *            the w weak
	 * @param is_valid
	 *            the is validx
	 */
	Correction(float west, float rest, float ww, float hh, float correction_, float wstrong,
			float wweak, bool is_ok) {
		w_est = west;
		r_est = rest;
		w = ww;
		h = hh;
		correction = correction_;
		w_strong = wstrong;
		w_weak = wweak;
		is_valid = is_ok;
	}
}
