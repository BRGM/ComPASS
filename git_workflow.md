git workflow for ComPASS development
====================================

:warning: *en français le temps de se mettre d'accord*

:warning: L'idée ici n'est pas de figer les choses mais plutôt de définir un socle de pratiques communes que l'on espère bonnes.

Ce qui suit est un mélange du retour d'expérience actuel sur charms, des conseils de Julien Bigot et Thomas Dufaud et d'informations glanés sur internet.

# Rappels

Il existe aujourd'hui deux projets sur le gitlab INRIA pour ComPASS :
- un projet *interne* correspondant à l'équipe ANR et rattaché au groupe
[charms](https://gitlab.inria.fr/charms) : on doit être identifié pour accéder au projet, on s'autorise le français pour discuter entre nous,
- un projet *public*, donc visible par le monde entier sans identification, rattaché au groupe [ComPASS](https://gitlab.inria.fr/compass) qui a vocation à être une vitrine externe : l'anglais est donc de rigueur et les commentaires/posts/issues... sont visibles de tous.

Il existe deux branches protégées sur le projet interne qui sont *master* et *develop*. La branche *master* suit la branche *develop* au fur et à mesure des releases. Elle a vocation à être utilisée en production et à aller sur le projet *public*. On vise a priori dans un premier temps des cycles de *release* assez fréquents pour que *master* et *develop* soient synchrones.

# Nature des interventions sur le code

On distinge trois type d'intervation, suivant leur durée:
- des modifications de type *typo* qui peuvent tout de même corriger des bugs importants
- des développements circonsrit autour d'une fonctionnalité (*feature*) 
- des développements de long terme *longterm* : refonte d'une partie du code, utilisation dans un cadre spécifique et mise au point de nouveaux cas test pour une thématique particulière

:question: Doit-on faire des commit sur *develop* pour des *typo* qui ne correspondent pas à la correction d'un bug

# Objectifs et régles de bases

On cherche à avoir un historique de développement le plus clair possible pour la branche *develop* dans un cadre collaboratif où plusieurs personnes peuvent intervenir en même temps sur les mêmes portions du code.

## Règles de base

Par importance décroissante:

- Se synchroniser le plus fréquement possible avec la branche develop.
- Faire des commits unitaires  bien documentés (cf. [section dédiée](#Faire-de-bons-commits)) : c'est particulièrement important pour pouvoir facilement identifier des régressions ou des bugs introduits grâce à `git bisect`.
- Mettre au propre son travail régulièrement (ou essayer de maintenir un écart raisonnable entre théorie et pratique sur ce point).
- De manière générale, demander de l'aide et communiquer sur son travail (en particulier en utilisant les outils *issue* et *milestone* de gitlab) même si c'est sans réponse :cry:, cela permet de garder une trace pour soi et pour les autres.

:question: Je rajouterais bien comme corrollaire à la première règle "et intégrer le plus souvent ses développements dans *develop* quitte à utiliser des if ou des blocs #ifdef" (cf. ce qui a été fait avec Laurence et l'option 
DEF_FREEFLOW_STRUCTURES dans le CMakeLists.txt).

# Mise en oeuvre

## Faire de bons commits

... ou en tout cas essayer !

Un commit doit être unitaire et doit rester centrer autour de l'introduction d'une fonctionnalité ou de la correction d'un bug. Diviser votre travail en interventions unitaires sur le code (ca vous aidera, si, si...) et faites autant de commit qu'il le faut.

Le fichier CodingConventions.md reprend des conseils standards pour effectuer de bon message de commit.

En particulier :
    * séparer sujet et corps du texte par une ligne blanche (en ligne de commande `git commit -m "Sujet" -m 'Corps du texte éventuellement sur plusieurs lignes'` produire le résultat escompté),
    * le sujet :
        + doit être clair, concis et précis, en anglais à la forme impérative
        + doit rester court, idéalement moins de 50 caractères, impérativement moins de 72
        + doit commencer par une majuscule et ne pas comporter de point final
    * le corps de texte :
        + ne doit pas comporter de lignes de plus de 72 caractères
        + doit expliquer ce que le commit fait et pourquoi, pas comment il le fait (si le commit est unitaire les lignes de code sont là pour expliquer le comment)
        + ne comporte pas d'indications temporelles (pensez à régler votre horloge système...)

Des balises en tête de sujet peuvent être utilisées pour identifier rapidement la nature du commit. Parmi celles déjà utilisées:
    - `BugFix` pour la correction d'un bug,
    - `WIP` (=*Work In Progress*) pour identifier un travail en cours dont la fusion sur la branche *develop* sera ainsi rendue impossible

Dans le message vous pouvez faire référence à :
  - une *issue*: #123
  - un *merge request*: !123
  - un *code snippet*: $123 (fonctionnalité pratique pour discuter de l'implémenation très spécifique).


Le commentaire du commit n'est pas nécessairement proportionnel au nombre de lignes de code modifiées.
Par un exemple une modification *typo* corrigeant un bug critique doit être documentée (et, si adapté, accompagnée d'un cas test qui aurait permis d'identifier le bug dans la pipeline d'intégration continue).

Des outils graphiques comme SourceTree permettent de faire facilement une sélection assez fine des modifications à introduire dans un commit (selection de lignes spécifiques, nettoyage des lignes inutiles...).

En résumé, que vous mettiez votre travail au propre ou que vous soyze dans une phase créatrice débridée, faire de *bons commits* vous simplifiera la vie. Si vous travaillez en local (ce qui est préférable) l'option `git --amend` permet de modifier un commit non poussé (cf. ci dessous mettre son travail au propre). 

:question: Doit-on préfixer le sujet des commits par des balises permettant de répérer quels parties du code ils concernent : je préférerais rattacher les commit à des issues auxquelles on peut affecter des étiquettes. La liste des fichiers affectés par le commit permet de voir rapidement quelle portion du code ils affectent.

## Mettre son travail au propre

Le contexte est celui ou on veut reconstruire un historique propre sur une branche *propre* dont le destion est d'être poussée sur le dépot pour un merge request. Cette historique est construit à partir de plusieurs branches *brouillons*.

Comme la phase de mise au propre d'un travail dépend de la méthode de travail retenue par chacun, il est difficile de donner des consignes générales mais quelques commandes git peuvent aider.

Dans ce qui suit __on suppose que l'on s'est mis sur la branche *propre*__ (i.e. `git checkout branche_propre`):

+ `git commit --amend` permet de modifier le dernier commit : travailler en local et **l'utiliser sans modération** que ce soit pour un commit sur une branche *brouillon* ou pour fignoler un commit sur la branche *propre* [doc](https://git-scm.com/docs/git-commit#Documentation/git-commit.txt---amend)
+ pour récupérer des modifications d'autres branches :
    ++ `git merge --no-commit` permet de fusionner une branche mais de s'arrêter avant de créer le commit et de pouvoir ainsi faire de nouvelles modifications [doc](https://git-scm.com/docs/git-merge#Documentation/git-merge.txt---no-commit)
    ++ `git merge --squash` fait la même chose que la commande précédent mais condense la fusion en un seul commit qui viendra s'insérer de manière linéaire dans l'historique de la branche propre [doc](https://git-scm.com/docs/git-merge#Documentation/git-merge.txt---squash)
    + `git cherry-pick -n` : permet de récupérer les modifications correspondant à un ou plusieurs commit(s) spécifique(s) mais sans créer de commit, qui permettra ensuite avec un `git commit` de créer un nouveau *sha-1* et de ne pas créer d'historiques de branches incohérents [doc](https://git-scm.com/docs/git-cherry-pick#Documentation/git-cherry-pick.txt--n).


Enfin pour les téméraires on peut utiliser `git rebase`, c'est pratique mais il faut surtout respecter la règle de **ne jamais faire un *rebase* sur une branche déjà poussée**.

> Quelques réflexions glanées sur internet à propos des différences merge/rebase :
>
> * Quand je fusionne une branche…
>     + Si elle est purement locale et temporaire, je m’assure qu’elle n’apparaît pas dans le graphe final de l’historique en faisant un fast-forward merge, ce qui peut nécessiter un rebase au préalable.
>     +  Si elle a une sémantique claire et documentée, je m’assure qu’elle apparaîtra clairement dans le graphe de l’historique, du début à la fin, en garantissant un true merge.
> * Quand je m’apprête à pusher mon travail local, je nettoie mon historique local d’abord pour partager un historique propre, au cordeau.
> * Quand je me vois refuser le push parce qu’un travail complémentaire a été pushé entre-temps, je rebase sur la branche distante à jour pour éviter de polluer le graphe par des tas de micro-merges malvenus.
>
> [source](https://delicious-insights.com/fr/articles/bien-utiliser-git-merge-et-rebase/) 

## Synchronisation avec *develop*

### *Merge fast-forward*

On a fait le choix de contraindre les merge vers develop à être *fast-forward* c'est à dire qu'il n'y a pas de commit correspondant au merge. Les modifications apportées sur la branche fusionnée doivent donc apparaitre comme ayant été faite après les dernières modifications de develop. La figure ci-dessous tirée de ce
[blog](http://www.dynamic-mess.com/developpement/git-merge-fast-forward-rebase-culture-git/) illustre cette différence avec à gauche la version *fast forward* versus la version avec un commit de merge à droite.

![différence fast forward no fast forward](http://www.dynamic-mess.com/Media/Images/git-ff.jpg "différence fast forward no fast forward")

> Concrètement, il s'agit d'intégrer l'historique de votre branche dans celui de la branche develop. Si vous aviez opté pour un merge sans fast-forward, vous auriez un historique distinct. Vous remarquerez qu'en l'absence de fast-forward, il y a un commit de merge (vert). J'aime bien cette approche, qui permet de conserver un historique complet et d'annuler le merge en seul coup ([source](http://www.dynamic-mess.com/developpement/git-merge-fast-forward-rebase-culture-git/))

Ce choix *force les développeurs à se synchroniser avec develop* et donc il les contraint à suivre une des règle de base. Il permet également de présenter historique linéaire de la branche *develop* : ce qui n'est par forcément un avantage si on ne veut pas avoir trop de commit sur develop. 

:exclamation: En particulier, avec cette contrainte on est sûr que si la pipeline d'intégration passe pour la branche fusionnée, elle passera aussi sur *develop*, en théorie, *_si_* la couverture par les tests est bonne on ne peut pas casser develop. Réciproquement si on casse develop c'est qu'il manque un test.

### En pratique

:exclamation::exclamation::exclamation: Tout sera toujours d'autant plus facile que les synchronisations avec *develop* sont fréquentes.

Pas mal de contraintes viennent du choix de n'avancer que par *fast-forward* sur develop sans avoir des messages de type "Merge feature 1 in branch develop".

:question: C'est peut être une complication inutile ? En particulier une stratégie qui peut marcher est de synchroniser sa branche avec develop en faisant un merge de *develop* dans sa branche mais cela créera dans l'historique de *develop* un messahe de type "Merge develop into vielle_branche" que l'on pourrait bien sur remplacer par "Synchronisation vielle_branche et develop".

#### Scénario 1 : intégration d'une nouvelle fonctionnalité

Considérons le scénario suivant :

> j'ai mis au propre le développement d'une nouvelle fonctionnalité sur *vieille_branche* et pendant le temps du développement et de la mise au propre je me suis desynchronisé de *develop*

Plusieurs options s'offrent à vous :
  
  + **vos modifications peuvent être condensées en un seul commit** :
    - vous vous synchronisez : `git pull develop`
    - vous vous mettez sur *vielle_branche* `git checkout vieille_branche`
    - vous faites un merge de develop : `git merge develop` (éventuellement `--no-commit` si vous voulez éditer le message et faire un commit manuel - ce qui vous sera proposé en cas de conflit)
        - cette étape va rajouter un commit de type "Merge de develop dans vielle_branche" que l'on ne veut pas avoir dans l'historique de develop
    - si tout semble bon en local vous poussez le nouvel état de *vieille_branche* sur le dépot
    - vous faites un merge request
        - un *maintainer* condensera votre commit en un seul commit en préservant la paternité de vos travaux et le fusionnera - il faudra certainenment interagir avec lui pour rédiger le message de commit qui va se retrouver dans develop (on peut utiliser pour cela l'interface du merge request de gitlab)
  
  + **vous avez une succession de commits que vous voulez préservez** :
    1 si votre perte de synchronisation avec develop est recente vous pouvez tenter un rebase automatique via l'interface de gitlab
        - si vous vous voulez être prudent créez une copie de vieill_branche (`git checkout -b vieill_branche_copy` en étant déjà sur *vieille_branche*)
        - vous poussez *vieille_branche* sur develop, faites un merge request et appuyez sur le bouton vert rebase [cf. la doc de gitlab](https://docs.gitlab.com/ee/user/project/merge_requests/fast_forward_merge.html#enabling-fast-forward-merges), si tout se passe bien et les tests d'integration passent c'est fini
    2 le *rebase* automatique n'aboutit pas
        - vous repartez en local de develop après vous être synchronisé `git pull develop`
        - plusieurs options pour recreer à la main un historique propre
            - soit vous tentez un git rebase interactif
            - soit vous créez à partir de develop une nouvelle branche `git checkout -b vielle_branche_merge` puis en utilisant les commandes git vous [recrivez un historique](#Mettre-son-travail-au-propre).

## Convention de nommage des branches

Je vous propose d'adopter la convention de nommage des branches suivante :
- pour les modifications de type *typos* un nom explicite ou la même convention de nommage que pour les features
- pour les modifications de type *feature*, créer d'abord un issue sur gitlab, même s'il n'y a qu'un titre préfixer ensuite toute les branches qui se rapportent à cette issue par son numéro, par exemple `123-what_and_why` pour une branche en lien avec l'issue 123. Le nom de la branche sera court mais idélament explicite. On peut utiliser les noms des développeurs lorsque l'on travail à plusieurs sur le même point
- pour les modifications de type *longterm*, créer d'abord un milestone dans gitlab et utiliser la balise `LT` et le numéro du milestone pour préfixer le nom de la branche, par exemple `LT123-revolution`. Lorsqu'une nouvelle fonctionnalité est insérée dans cette branche, suivre la convention précédente en rattachant simplement l'issue au milestone concerné dans gitlab, lors de la création de l'issue. L'utilisation de la balise `LT`permet d'éviter de supprimer une branche lors d'une opération de nettoyage.

Les branches *LT* doivent être considérées avec autant de soin que la branche *develop* et il vaut mieux les synchroniser régulièrement à develop. En particulier elles ne sont pas des branches de travail, mais des branches dans lesquelle on met au propre sont travail.

:question: Doit-on préfixer le nom des branches par des balises permettant de répérer quels parties du code elles concernent (c'est un peu le rôle des étiquettes que l'on peut associer au issue dans gitlab).

# Ce qui manque encore

pre-commit hook pour le formatage (cf. issue #1).

La mise en place de ces outils correspondra à un commit important modificant de nombreux fichiers.
C'est donc une opération à faire de manière concertée avec tous les développeurs intervenant sur le code, pour qu'ils repartent de la 
