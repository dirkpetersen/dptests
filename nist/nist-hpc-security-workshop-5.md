# 5th NIST High-Performance Computing Security Workshop Summary

## Executive Overview

The 5th NIST High-Performance Computing Security Workshop, held May 7-8, 2025, in Gaithersburg, MD, brought together government agencies, national laboratories, academic institutions, and industry partners to address the evolving cybersecurity challenges in HPC environments. The workshop focused on practical implementation of security frameworks, emerging threats, and collaborative approaches to securing critical research infrastructure.

## Key Themes and Developments

### NIST Standards Evolution

A major highlight was the publication of **NIST SP 800-234 (Initial Public Draft)** - the HPC Security Overlay - on May 1, 2025. This collaborative effort across 11 institutions provides tailored security controls built upon the moderate baseline defined in SP 800-53B. The overlay addresses 60 controls with HPC-specific guidance covering areas such as:

- Role-based access control in multi-user environments
- HPC-specific logging and monitoring
- User-developed software management
- Performance impact considerations
- Shared GPU and accelerator security

The **NIST Protecting CUI Series** also saw significant updates, with SP 800-171 Rev 3 and SP 800-171A Rev 3 published in May 2024, introducing organization-defined parameters and streamlined requirements for protecting Controlled Unclassified Information in nonfederal systems.

### HPC Security Technical Exchange (STX) Program

The HPC STX program has evolved from informal meetups between national laboratories to a structured government-wide community of interest. The 2024 STX brought together 80 registrants across government, contractors, and academia, covering topics from compliance baselines to incident handling. However, the planned 2025 STX (originally scheduled for April 1-4 at LLNL) was postponed to September 16-19 due to government travel restrictions that reduced attendance by half.

### Emerging Security Challenges

#### Container Security in HPC
Research presented by Dr. Yuede Ji revealed alarming vulnerabilities in HPC container images. A comprehensive study of 4,784 HPC container images found an average of **2,109 vulnerabilities per image**, with 65% classified as medium severity or higher. While 97% had low exploitation probability scores, the sheer volume of vulnerabilities presents significant risk. The research highlighted the need for:

- Better vulnerability scanners tailored for HPC environments
- Cross-language code similarity detection for vulnerability identification
- Secure container runtimes with hardware-based isolation

#### AI and Machine Learning Security
The integration of AI into scientific workflows introduces new attack vectors including:

- **Data poisoning** through tampered training datasets
- **Model drift** where validated models become less accurate over time
- **Agentic AI risks** with autonomous agents capable of lateral movement
- **Supply chain vulnerabilities** in shared models and datasets

Traditional cybersecurity approaches are insufficient for AI-enabled research environments, requiring new frameworks for securing data, models, and AI pipelines.

#### Kernel-Level Vulnerabilities
ORNL's research into HPC kernel fuzzing using Syzkaller demonstrated the productivity of targeted vulnerability research. Over two years of testing on the Cray Operating System led to 8 consistent crashes and 4 reproducible findings, with one vendor-accepted vulnerability assigned CVE-2025-27087. Scaling to virtual machine-based testing yielded 41 crashes and 1 reproducer in just three months.

### Compliance and Risk Management

#### Zero Trust Architecture
Multiple presentations addressed implementing zero trust principles in HPC environments. AWS outlined three building blocks for zero trust HPC:

1. **Hardware Root of Trust** through the Nitro System
2. **Identity & Access Management** with fine-grained policy controls
3. **VPC Micro Perimeters** providing dynamic isolation

The challenge remains balancing security granularity with operational complexity - determining whether isolation should occur at the cluster, queue, job, or task level.

#### Trusted CI Framework Adoption
The Trusted CI Framework has gained significant traction, with 17 organizations completing Framework Cohorts across five cohorts. The framework focuses on cybersecurity programmatics beyond technical controls, addressing mission alignment, governance, resources, and controls. ACCESS (NSF's cyberinfrastructure coordination hub) has adopted the framework to structure its cybersecurity program.

### Integrated Research Infrastructure (IRI)

DOE's vision for integrated research infrastructure represents a paradigm shift toward seamless, secure integration of experimental facilities, computing resources, and data repositories. The IRI approach addresses several architectural challenges:

- **Federated identity management** allowing single sign-on across facilities
- **Risk-based access controls** with different security tiers based on data sensitivity
- **Workflow orchestration** spanning diverse, distributed resources

Case studies like CITADEL at ORNL demonstrate successful implementation of secure computing environments for protected health information, achieving Authority to Operate (ATO) status for processing sensitive data on leadership-class systems.

## Implementation Challenges and Solutions

### Multi-Institutional Coordination
Successful HPC security implementation requires cross-institutional teams including:

- Office of Research
- Central and distributed IT
- Information security and privacy offices
- Legal counsel

Small institutions benefit from partnerships with larger "big brother/sister" institutions in their region for shared expertise and resources.

### Software Approval and User-Developed Code
The tension between security requirements and research flexibility remains a persistent challenge. NIST SP 800-234 acknowledges that "users may be allowed to install and develop software that is necessary for their mission," but implementing approved software repositories while maintaining research agility requires careful balance.

### Performance vs. Security Trade-offs
HPC environments face unique constraints where security controls that impact performance directly affect system cost and research productivity. Solutions must be designed with minimal performance overhead - AWS targets 2% impact for network security, while TEE-IO implementations typically see less than 5% performance degradation.

## Future Directions

### Supply Chain Security
Hardware and software supply chain risks are increasingly recognized as critical vulnerabilities. Microsoft's presentation on cryptographic approaches highlighted the need for:

- Device ownership transfer protocols
- Multiple endorsement chains (vendor + platform provider)
- Field entropy provisioning to protect against manufacturing compromises

### Quantum-Resistant Cryptography
The transition to quantum-resistant encryption algorithms presents both opportunities and challenges for HPC environments. While symmetric encryption impacts are expected to be minimal, asymmetric key exchanges may affect performance, and the migration timeline extends into 2027.

### Community Building and Knowledge Sharing
The workshop emphasized the importance of sustained community engagement through:

- Regular technical exchanges and working groups
- Shared vulnerability intelligence and threat information
- Collaborative development of security tools and frameworks
- Cross-sector partnerships between government, academia, and industry

## Conclusion

The 5th NIST HPC Security Workshop demonstrated the maturation of HPC cybersecurity as a distinct discipline requiring specialized approaches, standards, and community coordination. While significant progress has been made in developing frameworks like NIST SP 800-234 and building communities of practice, emerging challenges from AI integration, supply chain risks, and increasingly sophisticated threats require continued innovation and collaboration.

The postponement of the 2025 STX highlights the ongoing challenges of maintaining momentum in government-led initiatives, while the success of programs like the Trusted CI Framework Cohorts demonstrates the value of structured, collaborative approaches to cybersecurity program development. Moving forward, the community must balance the need for standardization with the flexibility required for cutting-edge research, ensuring that security enhances rather than impedes scientific discovery.

