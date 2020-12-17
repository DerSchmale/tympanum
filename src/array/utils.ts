export function removeElementOutOfOrder<T>(target: Array<T>, elm: T): number
{
    const last = target.pop();

    if (last === elm) {
        return target.length;
    }
    else {
        let index = target.indexOf(elm);

        if (index === -1) {
            target.push(last);
            throw new Error("Removing component that's not present");
        }

        target[index] = last;
        return index;
    }
}