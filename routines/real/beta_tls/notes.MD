# notes

## Kerr Related 
* Currently with no constraint on distance: 
    ```
    >>> len(np.where(df['Prob']>0.5)[0])
    12720
    >>> len(np.where(df['Prob']>0.75)[0])
    9267
    >>> len(np.where(df['Prob']>0.90)[0])
    7280
    >>> len(np.where(df['Prob']>0.99)[0])
    4728
    >>> probs = np.sort(df['Prob'])[::-1]
    >>> print(probs[9998])
    0.69
    >>> print(len(np.where(probs>0.69)[0]))
    9991
    >>> print(len(np.where(probs>0.68)[0]))
    10101
    ```

* Approximately the 68th percentile? 