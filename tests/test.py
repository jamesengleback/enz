import enz
import fire

def protein():
    prot = enz.Protein('4KEY.pdb', 
                        seq='MTIKEMPQPKTFGELKNLPLLNTDKPVQALMKIADELGEIFKFEAPGRVTRYLSSQRLIKEACDESRFDKNLSQALKFVRDFFGDGLVTSWTHEKNWKKAHNILLPSFSQQAMKGYHAMMVDIAVQLVQKWERLNADEHIEVPEDMTRLTLDTIGLCGFNYRFNSFYRDQPHPFITSMVRALDEAMNKLQRANPDDPAYDENKRQFQEDIKVMNDLVDKIIADRKASSDDLLTHMLNGKDPETGEPLDDENIRYQIITFLIAGHETTSGLLSFALYFLVKNPHVLQKAAEEAARVLVDPVPSYKQVKQLKYVGMVLNEALRLWPTAPAFSLYAKEDTVLGGEYPLEKGDELMVLIPQLHRDKTIWGDDVEEFRPERFENPSAIPQHAFKPFGNGQRACIGQQFALHEATLVLGMMLKHFDFEDHTNYELDIKETLTLKPEGFVVKAKSKKIPL',
                       keep=['HEM', 'not'])
    print(prot)
    return prot

def refold():
    prot = protein()
    prot.mutate(33,'A')
    prot.refold()
    return prot

def dock():
    prot = protein()
    results = prot.dock('CCCCCCCCCCCCCCCC=O', 
                        target_sites=[49,82,330,400,263,188],
                        exhaustiveness=1)
    print(results)
    print(len(results))
    print(results[0])
    print(results[0]['mol'].df)
    return results


def main():
    pass

if __name__ == '__main__':
    fire.Fire()
