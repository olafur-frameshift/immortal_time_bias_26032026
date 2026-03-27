> **Attach:** `figures/fig1_nmsc.png`

---

A 2013 study of 4.4 million Danes found that skin cancer patients had half the all-cause mortality of the general population. Huge dataset, near-complete follow-up, biologically plausible mechanism: more sun, more vitamin D, stronger bones, healthier cardiovascular system.

A randomly awarded lottery prize produces the same result.

The figure above shows a simulation of the study design. With no true causal effect built in, skin cancer patients show HR = 0.68 for all-cause mortality, closely mirroring the paper's reported HR = 0.52. The grey band is the culprit: immortal time.

Non-melanoma skin cancer is diagnosed in the 60s and 70s, well after the study's age-40 entry point. To be counted as a skin cancer patient, you had to survive to your diagnosis. Decades of guaranteed survival are silently attributed to the exposed group, while the general population (including everyone who died before a diagnosis was ever possible) pays the price in the comparison.

Replace the diagnosis with a prize awarded at random between ages 70 and 80: HR = 0.57 for all-cause mortality, no biological mechanism required. Start follow-up at the prize date instead? HR = 0.97. No effect.

→ The effect is a property of the design, not the exposure
→ Any post-baseline event that requires surviving to receive will produce this pattern
→ The fix: start follow-up at the index event, not study entry

Full write-up and reproducible R simulation: [LINK_PLACEHOLDER]

#Biostatistics #Epidemiology #ImmortalTimeBias #ObservationalResearch #Frameshift
